% DRAW_DOMINANT_ROOT_LOCUS
%
% IEEE 39节点系统的独立主导极点根轨迹脚本，参数与局部函数均包含在本文件。
% 这是随加载参数变化的特征值轨迹，不是调用rlocus扫描传递函数增益。
% 模型采用固定电压幅值、无损有功网络和VSG相角/频率动态。
%
% 计算流程：
% 1. 根据节点和线路数据建立关联矩阵B、线路权重a、平衡功率向量varpi。
% 2. 对每个lambda求B*diag(a)*sin(B'*theta)=lambda*varpi的平衡点。
% 3. 对平衡点线性化，消去零注入代数节点，构造状态矩阵A(lambda)。
% 4. 用eig(A)求极点，删除一个旋转零模态，选实部最大的非旋转极点。
% 5. 横轴画Re(s)，纵轴画Im(s)，颜色表示lambda/lambda_4。
% 6. 主图叠加若干条非主导近轴分支，但局部放大只保留s_cr主导分支。
%
% 计算从lambda=0继续至第一条线路相角达到pi/2的容量边界；本算例
% 通过末端极点校验确认到达s=0。代码没有计算Re(s)>0的平衡分支，
% 因而这张图证明“到达虚轴”，不能单凭它声称已穿过虚轴。

clearvars;
close all;
clc;
% 清空工作区、关闭旧图窗并清理命令行，保证每次运行从同一初始状态开始。

% -------------------------------------------------------------------------
% 输出目录和数值计算参数
% -------------------------------------------------------------------------
% 这些参数控制牛顿迭代、连续追踪、边界搜索和图形精度。
% 图形和数据默认保存在当前.m文件所在文件夹；交互执行时退回当前目录。
outputFolder = fileparts(mfilename('fullpath'));
if isempty(outputFolder)
    outputFolder = pwd;
end

baseMVA = 100;                       % 功率标幺基准，单位MVA。
sourceTolerance = 1e-8;              % 判断节点净注入是否为零的容差。
newtonTolerance = 1e-11;             % 平衡方程无穷范数残差的收敛阈值。
newtonMaximumIterations = 150;       % 单次平衡点牛顿迭代上限。
lineSearchMaximumIterations = 60;    % 步长减半的最大次数。
continuationSubsteps = 20;           % 相邻加载点间继续追踪的子步数。
boundaryBisectionIterations = 50;    % 边界区间二分次数。
boundaryCandidateEdges = 8;          % 检验相角最大的前8条候选线路。
boundaryExpansionFactor = 2;       % 搜索上界时lambda的放大倍数。
boundaryMaximumLoadingFactor = 100;  % 边界搜索的加载上限。
boundaryAngle = pi/2;                % 相角容量边界，单位rad。
phaseAngleMargin = 1e-8;             % 普通牛顿求解需留在pi/2以内的裕量。
boundaryResidualTolerance = 5e-9;    % 增广边界方程残差容差。
boundaryAngleTolerance = 5e-7;       % 边界相角误差容差，单位rad。

% 论文中的参考量：lambda_4用于颜色条归一化，其他量用于构造VSG参数。
lambda4 = 4.15311915166592;             % 论文参考加载量，只用于归一化标注。
nominalFrequencyHz = 50;               % 额定频率记录值，未直接进入本状态矩阵。
frequencyDropHzAtLambda1 = 0.20;        % 用于确定总阻尼的频率偏移，单位Hz。
dampingShapeRange = [0.60 1.40];        % 相对阻尼分配范围，无量纲。
inertiaTimeConstantRange = [0.20 1.00]; % M_i/D_i的范围，单位s。

% -------------------------------------------------------------------------
% IEEE 39节点系统数据
% -------------------------------------------------------------------------
% 下面直接写入电压、负荷、发电机和线路参数，使脚本可以独立运行。
% 电压幅值和有功负荷来自case39；两列向量按节点1到39排列。
% 电压幅值单位为p.u.，在lambda扫描中固定不变。
busVoltageMagnitude = [ ...
    1.0393836; 1.0484941; 1.0307077; 1.0044600; 1.0060063; ...
    1.0082256; 0.99839728; 0.99787232; 1.0383320; 1.0178431; ...
    1.0133858; 1.0008150; 1.0149230; 1.0123190; 1.0161854; ...
    1.0325203; 1.0342365; 1.0315726; 1.0501068; 0.99101054; ...
    1.0323192; 1.0501427; 1.0451451; 1.0380010; 1.0576827; ...
    1.0525613; 1.0383449; 1.0503737; 1.0501149; 1.0499000; ...
    0.9820000; 0.9841000; 0.9972000; 1.0123000; 1.0494000; ...
    1.0636000; 1.0275000; 1.0265000; 1.0300000];

% 节点有功负荷单位MW，正值表示消耗。
busActiveDemandMW = [ ...
     97.6;   0.0; 322.0; 500.0;   0.0;   0.0; 233.8; 522.0; ...
      6.5;   0.0;   0.0;   8.53;  0.0;   0.0; 320.0; 329.0; ...
      0.0; 158.0;  0.0; 680.0; 274.0;   0.0; 247.5; 308.6; ...
    224.0; 139.0; 281.0; 206.0; 283.5;   0.0;   9.2;   0.0; ...
      0.0;   0.0;   0.0;   0.0;   0.0;   0.0; 1104.0];

% 在运发电机位于节点30至39，功率单位MW，与节点向量逐项对应。
generatorBus = (30:39).';
generatorActivePowerMW = [ ...
    250.000; 677.871; 650.000; 632.000; 508.000; ...
    650.000; 560.000; 540.000; 830.000; 1000.000];

% 支路各列依次为[起始节点,终止节点,电抗p.u.,变压器变比]。
% 忽略电阻和无功/电压动态；原MATPOWER数据的零变比按1处理。
branchData = [ ...
     1   2  0.0411  1.000;
     1  39  0.0250  1.000;
     2   3  0.0151  1.000;
     2  25  0.0086  1.000;
     2  30  0.0181  1.025;
     3   4  0.0213  1.000;
     3  18  0.0133  1.000;
     4   5  0.0128  1.000;
     4  14  0.0129  1.000;
     5   6  0.0026  1.000;
     5   8  0.0112  1.000;
     6   7  0.0092  1.000;
     6  11  0.0082  1.000;
     6  31  0.0250  1.070;
     7   8  0.0046  1.000;
     8   9  0.0363  1.000;
     9  39  0.0250  1.000;
    10  11  0.0043  1.000;
    10  13  0.0043  1.000;
    10  32  0.0200  1.070;
    12  11  0.0435  1.006;
    12  13  0.0435  1.006;
    13  14  0.0101  1.000;
    14  15  0.0217  1.000;
    15  16  0.0094  1.000;
    16  17  0.0089  1.000;
    16  19  0.0195  1.000;
    16  21  0.0135  1.000;
    16  24  0.0059  1.000;
    17  18  0.0082  1.000;
    17  27  0.0173  1.000;
    19  20  0.0138  1.060;
    19  33  0.0142  1.070;
    20  34  0.0180  1.009;
    21  22  0.0140  1.000;
    22  23  0.0096  1.000;
    22  35  0.0143  1.025;
    23  24  0.0350  1.000;
    23  36  0.0272  1.000;
    25  26  0.0323  1.000;
    25  37  0.0232  1.025;
    26  27  0.0147  1.000;
    26  28  0.0474  1.000;
    26  29  0.0625  1.000;
    28  29  0.0151  1.000;
    29  38  0.0156  1.025];

% 从数据维度读取节点数、线路数和各线路参数。
N = numel(busVoltageMagnitude);
numberOfEdges = size(branchData,1);
fromBus = branchData(:,1);
toBus = branchData(:,2);
reactance = branchData(:,3);
tapRatio = branchData(:,4);

% -------------------------------------------------------------------------
% 构造定电压无损网络和功率加载方向
% -------------------------------------------------------------------------
% incidence是节点-支路关联矩阵；lineWeight是线路传输系数。
% B的每列表示一条有向线路：起点+1、终点-1，故B'*theta为线路相角差。
incidence = sparse(N,numberOfEdges);
for edge = 1:numberOfEdges
    incidence(fromBus(edge),edge) = 1;
    incidence(toBus(edge),edge) = -1;
end

% a_e=V_i*V_j/(x_e*abs(tap_e))，支路有功为a_e*sin(theta_i-theta_j)。
lineWeight = busVoltageMagnitude(fromBus).*busVoltageMagnitude(toBus) ./ ...
    (reactance.*abs(tapRatio));

% 把发电功率累加到各节点并标幺化，节点净注入=发电功率-负荷功率。
generatorPowerPU = accumarray(generatorBus, ...
    generatorActivePowerMW/baseMVA,[N 1],@sum,0);
activeDemandPU = busActiveDemandMW/baseMVA;
netInjectionAC = generatorPowerPU-activeDemandPU;

% 正注入节点和负注入节点分别统计供给/需求总量，确定无损平衡缩放比例。
positiveInjection = netInjectionAC > sourceTolerance;
negativeInjection = netInjectionAC < -sourceTolerance;
sourceTotalPU = sum(netInjectionAC(positiveInjection));
demandTotalPU = -sum(netInjectionAC(negativeInjection));
sourceRescaling = demandTotalPU/sourceTotalPU;

% 原始交流系统含有损耗，先按比例缩放正注入，使无损模型满足功率平衡。
varpiNominal = netInjectionAC;
varpiNominal(positiveInjection) = ...
    sourceRescaling*netInjectionAC(positiveInjection);
varpiNominal(abs(varpiNominal) <= sourceTolerance) = 0;

% 动态模型把所有非零注入节点作为VSG节点，零注入节点作为代数节点。
vsgBuses = find(abs(varpiNominal) > sourceTolerance);
algebraicBuses = find(abs(varpiNominal) <= sourceTolerance);
referenceBus = 31;  % 固定theta_31=0，只消除平衡求解中的统一转角自由度。

% 根据节点功率份额构造非均匀阻尼和惯量；功率份额越大，参数位置越靠近对应范围端点。
absoluteShare = abs(varpiNominal(vsgBuses));
absoluteShare = absoluteShare/sum(absoluteShare);
shareRange = max(absoluteShare)-min(absoluteShare);
% 将功率份额映射到[0,1]；份额相同时取中点，避免除零。
if shareRange < eps
    ratingPosition = 0.5*ones(size(absoluteShare));
else
    ratingPosition = (absoluteShare-min(absoluteShare))/shareRange;
end

% D总量=总需求/角频率偏移，再按dampingShape分配给各VSG。
frequencyDropRadPerSecond = 2*pi*frequencyDropHzAtLambda1;
dampingTotal = demandTotalPU/frequencyDropRadPerSecond;
dampingShape = dampingShapeRange(1) + ...
    diff(dampingShapeRange)*ratingPosition;
damping = dampingTotal*dampingShape/sum(dampingShape);

% 先分配M_i/D_i时间常数，再通过M_i=D_i*(M_i/D_i)计算惯量。
inertiaTimeConstant = inertiaTimeConstantRange(1) + ...
    diff(inertiaTimeConstantRange)*(1-ratingPosition);
inertia = damping.*inertiaTimeConstant;

% 在继续求解平衡支路前，检查网络连通性和加载向量的零和条件。
assert(rank(full(incidence),1e-10) == N-1, ...
    'The embedded IEEE 39-bus network is disconnected.');
assert(abs(sum(varpiNominal)) < 1e-10, ...
    'The embedded balanced loading vector is not zero-sum.');

% -------------------------------------------------------------------------
% 搜索第一条相角容量边界
% -------------------------------------------------------------------------
% 先扩展lambda寻找可行/不可行区间，再用增广方程精确求出
% 某条线路相角达到pi/2时的lambdaCrossing和thetaCrossing。
% 把求解器参数装入结构体，传给本文件末尾的局部求解函数。
boundarySettings.newtonTolerance = newtonTolerance;
boundarySettings.newtonMaximumIterations = newtonMaximumIterations;
boundarySettings.lineSearchMaximumIterations = lineSearchMaximumIterations;
boundarySettings.continuationSubsteps = continuationSubsteps;
boundarySettings.boundaryBisectionIterations = boundaryBisectionIterations;
boundarySettings.boundaryCandidateEdges = boundaryCandidateEdges;
boundarySettings.boundaryExpansionFactor = boundaryExpansionFactor;
boundarySettings.boundaryMaximumLoadingFactor = ...
    boundaryMaximumLoadingFactor;
boundarySettings.boundaryAngle = boundaryAngle;
boundarySettings.phaseAngleMargin = phaseAngleMargin;
boundarySettings.boundaryResidualTolerance = boundaryResidualTolerance;
boundarySettings.boundaryAngleTolerance = boundaryAngleTolerance;

[lambdaCrossing,thetaCrossing,limitingEdge,limitingSign, ...
    boundaryResidual,boundaryAngleError] = find_phase_boundary_standalone( ...
    incidence,lineWeight,varpiNominal,referenceBus,boundarySettings);

% 在0到边界之间取密集加载点，并追加精确边界点。
% 边界点包含即将到达虚轴的非旋转零特征值；每个加载点只删除统一转角
% 产生的旋转零模态。临界点附近采用对数加密，避免s_cr局部图只剩一个点。
numberOfLoadingPoints = 241; % 远离临界点的基础采样点数。
preBoundaryGap = 1e-6;       % 最后一个普通点距临界lambda的最小距离。
boundarySamplingWidth = 0.02; % 主图保留的临界前完整过渡区间。
boundaryZoomWidth = 2e-4;     % 局部放大框只显示临界点附近的窄区间。
numberOfBoundaryZoomPoints = 160; % 临界前过渡区间的对数采样点数。
lambdaBase = linspace(0,max(0,lambdaCrossing-boundarySamplingWidth), ...
    numberOfLoadingPoints).';
% 主图从较早的主导实极点开始保留，逐渐加密到虚轴临界点。
boundaryGaps = logspace(log10(boundarySamplingWidth), ...
    log10(preBoundaryGap),numberOfBoundaryZoomPoints).';
lambdaNearBoundary = lambdaCrossing-boundaryGaps;
lambdaValues = unique([lambdaBase;lambdaNearBoundary;lambdaCrossing]);
lambdaValues = sort(lambdaValues);
numberOfLoadingPoints = numel(lambdaValues);

% NaN表示该加载点未选到这一类别的主导极点，绘图时不会跨缺失点连线。
numberOfVSG = numel(vsgBuses);
dominantPositive = NaN(numberOfLoadingPoints,1);
dominantNegative = NaN(numberOfLoadingPoints,1);
dominantReal = NaN(numberOfLoadingPoints,1);
maximumEdgeAngle = NaN(numberOfLoadingPoints,1);
powerFlowResidual = NaN(numberOfLoadingPoints,1);
maximumRealPart = NaN(numberOfLoadingPoints,1);
% 保存全部非旋转极点，并在相邻加载点之间做最近距离模态匹配。
% 这样绘制全部分支时不会因极点按实部排序变化而产生错误的竖直连线。
numberOfNonrotationalPoles = 2*numberOfVSG-1;
allNonrotationalPoles = NaN(numberOfNonrotationalPoles,numberOfLoadingPoints);
previousNonrotationalPoles = [];

% lambda=0时各节点相角相同，取全零作为初始平衡点。
theta = zeros(N,1);
lambdaPrevious = 0;
phaseLimit = pi/2-phaseAngleMargin;

for point = 1:numberOfLoadingPoints
    lambda = lambdaValues(point);
    isBoundaryPoint = point == numberOfLoadingPoints;

    % 逐点继续求解潮流平衡点，并在每个平衡点线性化计算系统极点。

    % 临界点直接使用增广方程解，其余点以上一点为初值进行继续追踪。
    if isBoundaryPoint
        theta = thetaCrossing;
    elseif point > 1
        [theta,solver] = continue_equilibrium_standalone( ...
            incidence,lineWeight,varpiNominal,referenceBus, ...
            lambdaPrevious,theta,lambda,boundarySettings, ...
            phaseLimit);
        if ~solver.converged
            error(['The standalone continuation failed at lambda = ', ...
                '%.12g: %s'],lambda,solver.reason);
        end
    end

    [poles,metrics] = reduced_vsg_poles_standalone( ...
        incidence,lineWeight,varpiNominal,theta,lambda, ...
        vsgBuses,algebraicBuses,damping,inertia);
    dominant = select_dominant_poles_standalone(poles);

    % 删除旋转零模态，并将全部其余极点按相邻加载点的最近距离匹配。
    polesForLocus = poles;
    [~,rotationIndex] = min(abs(polesForLocus));
    polesForLocus(rotationIndex) = [];
    if point == 1
        [~,poleOrder] = sortrows([real(polesForLocus),imag(polesForLocus)], ...
            [1 2]);
        orderedNonrotationalPoles = polesForLocus(poleOrder);
    else
        orderedNonrotationalPoles = match_poles_by_distance_standalone( ...
            previousNonrotationalPoles,polesForLocus);
    end
    allNonrotationalPoles(:,point) = orderedNonrotationalPoles;
    previousNonrotationalPoles = orderedNonrotationalPoles;

    % 保存本次选中的上复极点、下复极点或实极点，以及数值质量指标。
    dominantPositive(point) = dominant.positive;
    dominantNegative(point) = dominant.negative;
    dominantReal(point) = dominant.real;
    maximumEdgeAngle(point) = metrics.maximumEdgeAngle;
    powerFlowResidual(point) = metrics.powerFlowResidual;
    maximumRealPart(point) = dominant.maximumRealPart;
    lambdaPrevious = lambda;
end

% 检查最后一个数据点确实是虚轴临界点。
assert(abs(lambdaValues(end)-lambdaCrossing) < 1e-12);
assert(isfinite(dominantReal(end)) && ...
    abs(dominantReal(end)) < 1e-6, ...
    'The non-rotational pole did not reach the imaginary axis.');

% -------------------------------------------------------------------------
% 保存根轨迹数值数据
% -------------------------------------------------------------------------
% CSV便于查看，MAT文件保存后续分析所需的完整变量。
% 表中Re/Im为极点实部/虚部，单位s^-1；相角单位rad。
loadingRatio = lambdaValues/lambda4;
dominantTable = table(lambdaValues,loadingRatio, ...
    real(dominantPositive),imag(dominantPositive), ...
    real(dominantNegative),imag(dominantNegative), ...
    dominantReal,maximumRealPart,maximumEdgeAngle,powerFlowResidual, ...
    'VariableNames',{'LoadingFactor','LoadingRatio', ...
    'PositivePoleRealPart','PositivePoleImagPart', ...
    'NegativePoleRealPart','NegativePoleImagPart', ...
    'RealDominantPole','MaximumRealPart','MaximumEdgeAngle', ...
    'PowerFlowResidual'});
% CSV只保存各扫描点的数值结果，不是脚本运行所需的输入文件。
writetable(dominantTable,fullfile(outputFolder, ...
    'dominant_root_locus_data.csv'));

% MAT附带网络参数、临界线路和阻尼/惯量，便于后续复核。
save(fullfile(outputFolder,'dominant_root_locus_data.mat'), ...
    'dominantTable','lambda4','lambdaCrossing','limitingEdge', ...
    'limitingSign','boundaryResidual','boundaryAngleError', ...
    'vsgBuses','algebraicBuses','damping','inertia', ...
    'varpiNominal','lineWeight','branchData','busVoltageMagnitude', ...
    'allNonrotationalPoles');

% -------------------------------------------------------------------------
% 绘制论文风格的主导极点根轨迹
% -------------------------------------------------------------------------
% 主图显示上、下共轭极点和接近虚轴的实极点；颜色表示lambda/lambda_4。
% 只创建这一张活动图窗；Position分别为[左,下,宽,高]，单位inch。
figureHandle = figure('Color','w','Units','inches', ...
    'Position',[1 1 6.50 4.80]);
% 主坐标区位置采用整张图的归一化比例。
axesHandle = axes(figureHandle,'Position',[0.12 0.12 0.70 0.78]);
hold(axesHandle,'on');

colorMap = loading_colormap_standalone(256);
colormap(figureHandle,colorMap);
% 把色标上限向上取至0.2的整数倍，并在主图/各放大框中保持相同范围。
colorLimit = max(1,ceil(max(loadingRatio)*5)/5);

% 灰色细线连接轨迹，彩色散点编码加载比例，风格与参考图一致。
plot(axesHandle,real(dominantPositive),imag(dominantPositive), ...
    '-','Color',[0.25 0.25 0.25],'LineWidth',0.75, ...
    'HandleVisibility','off');
plot(axesHandle,real(dominantNegative),imag(dominantNegative), ...
    '-','Color',[0.25 0.25 0.25],'LineWidth',0.75, ...
    'HandleVisibility','off');
plot(axesHandle,real(dominantReal),imag(dominantReal), ...
    '-','Color',[0.25 0.25 0.25],'LineWidth',0.85, ...
    'HandleVisibility','off');

scatter(axesHandle,real(dominantPositive),imag(dominantPositive), ...
    18,loadingRatio,'filled','MarkerEdgeColor',[0.12 0.12 0.12], ...
    'LineWidth',0.35, ...
    'HandleVisibility','off');
scatter(axesHandle,real(dominantNegative),imag(dominantNegative), ...
    18,loadingRatio,'filled','MarkerEdgeColor',[0.12 0.12 0.12], ...
    'LineWidth',0.35, ...
    'HandleVisibility','off');
scatter(axesHandle,real(dominantReal),imag(dominantReal), ...
    22,loadingRatio,'filled','MarkerEdgeColor',[0.12 0.12 0.12], ...
    'LineWidth',0.35, ...
    'HandleVisibility','off');

% 用细线和小点显示全部非旋转极点分支；主导极点随后用大彩色点覆盖。
for branchIndex = 1:numberOfNonrotationalPoles
    branchX = real(allNonrotationalPoles(branchIndex,:));
    branchY = imag(allNonrotationalPoles(branchIndex,:));
    plot(axesHandle,branchX,branchY,'-', ...
        'Color',[0.68 0.68 0.68],'LineWidth',0.28, ...
        'HandleVisibility','off');
    scatter(axesHandle,branchX,branchY,6,loadingRatio, ...
        'filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.28, ...
        'HandleVisibility','off');
end

% 用空心圆标出非旋转主导极点到达虚轴的位置s=0。
plot(axesHandle,0,0,'ko','MarkerSize',4.5,'LineWidth',0.8, ...
    'MarkerFaceColor','w','HandleVisibility','off');
xline(axesHandle,0,'k-','LineWidth',0.85,'HandleVisibility','off');
yline(axesHandle,0,'k-','LineWidth',0.70,'HandleVisibility','off');

% 统一字体、边框、刻度方向和浅色虚线网格。
set(axesHandle,'FontName','Times New Roman','FontSize',10.5, ...
    'FontWeight','bold','LineWidth',0.7,'TickDir','out', ...
    'Box','on','Layer','top');
grid(axesHandle,'on');
axesHandle.GridAlpha = 0.16;
axesHandle.GridLineStyle = ':';

% 主图范围覆盖所有已选主导极点，并在虚轴右侧留少量空白。
allRealParts = [real(dominantPositive);real(dominantNegative); ...
    real(dominantReal);real(allNonrotationalPoles(:))];
allImagParts = [imag(dominantPositive);imag(dominantNegative); ...
    imag(dominantReal);imag(allNonrotationalPoles(:))];
allRealParts = allRealParts(isfinite(allRealParts));
allImagParts = allImagParts(isfinite(allImagParts));
xLower = min(allRealParts)-0.08;
xUpper = 0.05;
yLimit = max(1,max(abs(allImagParts)))*1.10;
xlim(axesHandle,[xLower xUpper]);
ylim(axesHandle,[-yLimit yLimit]);
pbaspect(axesHandle,[1.25 1 1]);

xlabel(axesHandle,'Real part (s^{-1})', ...
    'FontName','Times New Roman','FontSize',12,'FontWeight','bold');
ylabel(axesHandle,'Imaginary part (s^{-1})', ...
    'FontName','Times New Roman','FontSize',12,'FontWeight','bold');

% 独立放置色条及lambda/lambda_4标签，避免色条压缩局部框。
colorbarHandle = colorbar(axesHandle);
set(colorbarHandle,'Units','normalized','Position',[0.87 0.12 0.035 0.62], ...
    'FontName','Times New Roman','FontSize',10,'FontWeight','bold', ...
    'LineWidth',0.7);
clim(axesHandle,[0 colorLimit]);
colorbarHandle.Ticks = unique([0:0.2:colorLimit,colorLimit]);
annotation(figureHandle,'textbox',[0.84 0.79 0.12 0.07], ...
    'String','\lambda/\lambda_4','Interpreter','tex', ...
    'FontName','Times New Roman','FontSize',10,'FontWeight','bold', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'Rotation',0,'LineStyle','none','FitBoxToText','off');
% 用一个黑边彩色圆点说明主图中较大的黑边散点就是主导极点。
annotation(figureHandle,'ellipse',[0.66 0.910 0.012 0.020], ...
    'FaceColor',[0.85 0.30 0.10], ...
    'EdgeColor',[0.10 0.10 0.10],'LineWidth',0.6);
annotation(figureHandle,'textbox',[0.67 0.903 0.15 0.035], ...
    'String','Dominant poles','Interpreter','none', ...
    'FontName','Times New Roman','FontSize',9,'FontWeight','bold', ...
    'Color',[0.10 0.10 0.10], ...
    'HorizontalAlignment','left','VerticalAlignment','middle', ...
    'LineStyle','none','FitBoxToText','off');

% s_cr放大框只显示临界主导实极点；不再放大非临界复极点分支。
zoomXLimits = [-0.10 0.00];   % 局部实部范围，单位s^-1。
zoomYLimits = [-0.0010 0.0010];% 局部虚部范围，单位s^-1。
finiteRealPole = isfinite(dominantReal);
% 局部框只取临界点附近的主导实极点；主图仍保留完整过渡分支。
zoomPoleMask = finiteRealPole & ...
    (lambdaValues >= lambdaCrossing-boundaryZoomWidth);
zoomReal = dominantReal(zoomPoleMask);
zoomRatio = loadingRatio(zoomPoleMask);

pngFile = fullfile(outputFolder,'dominant_root_locus.png');
pdfFile = fullfile(outputFolder,'dominant_root_locus.pdf');
% The full figure is exported after the three inset axes have been added.

% 已注释：独立临界点放大图。保留旧代码供参考，不创建图窗也不导出。
%{
% Export the local view separately for direct use in the paper.
zoomFigure = figure('Color','w','Units','inches', ...
    'Position',[1 1 4.25 3.25]);
colormap(zoomFigure,colorMap);
zoomAxesStandalone = axes(zoomFigure,'Position',[0.14 0.18 0.70 0.68]);
hold(zoomAxesStandalone,'on');
yline(zoomAxesStandalone,0,'k-','LineWidth',0.60,'HandleVisibility','off');
xline(zoomAxesStandalone,0,'k-','LineWidth',0.75,'HandleVisibility','off');
plot(zoomAxesStandalone,zoomReal,zeros(size(zoomReal)),'-', ...
    'Color',[0.25 0.25 0.25],'LineWidth',1.0,'HandleVisibility','off');
scatter(zoomAxesStandalone,zoomReal,zeros(size(zoomReal)),28,zoomRatio, ...
    'filled','MarkerEdgeColor','none','HandleVisibility','off');
plot(zoomAxesStandalone,0,0,'ko','MarkerSize',5,'LineWidth',0.85, ...
    'MarkerFaceColor','w','HandleVisibility','off');
xlim(zoomAxesStandalone,zoomXLimits);
ylim(zoomAxesStandalone,zoomYLimits);
set(zoomAxesStandalone,'FontName','Times New Roman','FontSize',10.5, ...
    'LineWidth',0.7,'TickDir','out','Box','on','Layer','top');
grid(zoomAxesStandalone,'on');
zoomAxesStandalone.GridAlpha = 0.16;
zoomAxesStandalone.GridLineStyle = ':';
xlabel(zoomAxesStandalone,'Real part (s^{-1})', ...
    'FontName','Times New Roman','FontSize',10);
ylabel(zoomAxesStandalone,'Imaginary part (s^{-1})', ...
    'FontName','Times New Roman','FontSize',10);
title(zoomAxesStandalone,'Virtual-axis crossing', ...
    'FontName','Times New Roman','FontSize',10,'FontWeight','normal');
zoomColorbar = colorbar(zoomAxesStandalone);
set(zoomColorbar,'Units','normalized','Position',[0.87 0.27 0.035 0.50], ...
    'FontName','Times New Roman','FontSize',10,'LineWidth',0.7);
clim(zoomAxesStandalone,[0 colorLimit]);
zoomColorbar.Ticks = unique([0:0.2:colorLimit,colorLimit]);
annotation(zoomFigure,'textbox',[0.83 0.79 0.13 0.07], ...
    'String','$\lambda/\lambda_4$','Interpreter','latex', ...
    'FontName','Times New Roman','FontSize',10, ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'LineStyle','none','FitBoxToText','off');

zoomPngFile = fullfile(outputFolder,'dominant_root_locus_zoom.png');
zoomPdfFile = fullfile(outputFolder,'dominant_root_locus_zoom.pdf');
try
    exportgraphics(zoomFigure,zoomPngFile,'Resolution',300);
    exportgraphics(zoomFigure,zoomPdfFile,'ContentType','vector');
catch zoomExportException
    warning('The standalone crossing zoom could not be exported: %s', ...
        zoomExportException.message);
end
%}

% 提取有限的上、下复极点，并根据数据范围设置两个局部坐标窗口。
finitePositive = isfinite(dominantPositive);
finiteNegative = isfinite(dominantNegative);
complexRealParts = [real(dominantPositive(finitePositive)); ...
    real(dominantNegative(finiteNegative))];
complexImagParts = [imag(dominantPositive(finitePositive)); ...
    imag(dominantNegative(finiteNegative))];
complexXLimits = [min(complexRealParts)-0.002, ...
    max(complexRealParts)+0.002];
positiveImagLimits = [min(imag(dominantPositive(finitePositive)))-0.03, ...
    max(imag(dominantPositive(finitePositive)))+0.03];
negativeImagLimits = [min(imag(dominantNegative(finiteNegative)))-0.03, ...
    max(imag(dominantNegative(finiteNegative)))+0.03];

% 在同一张主图中加入临界主导实极点的局部放大框。
% 非临界的s_+、s_-局部放大框暂时注释掉，只保留s_cr。
insetPositions = [ ...
    0.18 0.61 0.25 0.22];      % 临界点放大框，向右移开主图中的-0.5刻度。
% 每个放大框都在同一figure中建立独立axes，并复用同一批计算数据。
insetAxes = gobjects(1,1);
for insetIndex = 1
    insetAxes(insetIndex) = axes(figureHandle, ...
        'Position',insetPositions(insetIndex,:), ...
        'Color','w','Box','on');
    hold(insetAxes(insetIndex),'on');
    % 这里只放大临界主导实极点s_cr。
    insetX = zoomReal;
    insetY = zeros(size(zoomReal));
    insetRatio = zoomRatio;
    insetXLimits = zoomXLimits;
    insetYLimits = zoomYLimits;
    insetTitle = '\bf s_{\rm cr}';

    % 只有近零实极点框需要画x=0、y=0两条轴线。
    if insetIndex == 1
        yline(insetAxes(insetIndex),0,'k-','LineWidth',0.55, ...
            'HandleVisibility','off');
        xline(insetAxes(insetIndex),0,'k-','LineWidth',0.65, ...
            'HandleVisibility','off');
    end
    plot(insetAxes(insetIndex),insetX,insetY,'-', ...
        'Color',[0.25 0.25 0.25],'LineWidth',0.8,'HandleVisibility','off');
    scatter(insetAxes(insetIndex),insetX,insetY,20,insetRatio, ...
        'filled','MarkerEdgeColor',[0.12 0.12 0.12], ...
        'LineWidth',0.30,'HandleVisibility','off');
    if insetIndex == 1
        plot(insetAxes(insetIndex),0,0,'ko','MarkerSize',3.8, ...
            'LineWidth',0.65,'MarkerFaceColor','w','HandleVisibility','off');
    end
    % 每个框的x/y范围独立，但颜色范围与主图完全相同。
    xlim(insetAxes(insetIndex),insetXLimits);
    ylim(insetAxes(insetIndex),insetYLimits);
    set(insetAxes(insetIndex),'FontName','Times New Roman','FontSize',8, ...
        'FontWeight','bold','LineWidth',0.5,'TickDir','out', ...
        'Layer','top');
    grid(insetAxes(insetIndex),'on');
    clim(insetAxes(insetIndex),[0 colorLimit]);
    insetAxes(insetIndex).GridAlpha = 0.16;
    insetAxes(insetIndex).GridLineStyle = ':';
    % 保留临界主导分支的刻度和坐标标签。
    % xlabel(insetAxes(insetIndex),'Re(s)', ...
    %     'FontName','Times New Roman','FontSize',8);
    % ylabel(insetAxes(insetIndex),'Im(s)', ...
    %     'FontName','Times New Roman','FontSize',8);
    % title(insetAxes(insetIndex),insetTitle,'FontName','Times New Roman', ...
    %     'FontSize',11,'FontWeight','bold','Interpreter','latex');
    % 只在临界点框标注绝对lambda及相对lambda_4比例。
    if insetIndex == 1
        text(insetAxes(insetIndex),insetXLimits(1)+0.001, ...
            insetYLimits(2)-0.00005, ...
            ['\lambda_{cr}=',num2str(lambdaCrossing,'%.6f')], ...
            'FontName','Times New Roman','FontSize',10,'FontWeight','bold', ...
            'VerticalAlignment','top','Interpreter','tex');
        text(insetAxes(insetIndex),insetXLimits(1)+0.001, ...
            insetYLimits(2)-0.00045, ...
            ['\lambda_{cr}/\lambda_4=', ...
            num2str(lambdaCrossing/lambda4,'%.6f')], ...
            'FontName','Times New Roman','FontSize',10,'FontWeight','bold', ...
            'VerticalAlignment','top','Interpreter','tex');
    end
end

% % ========================================================================
% % 临界主导分支箭头方向和位置设置
% % ========================================================================
% % annotation('arrow',[x_start,x_end],[y_start,y_end])的箭头尖端在
% % [x_end,y_end]，箭尾在[x_start,y_start]。因此：
% %   1) 想反向：交换arrowSources和arrowTargets的两个坐标；
% %   2) 想移动箭尾：修改crossingSource；
% %   3) 想移动箭头尖端：修改crossingTarget；
% %   4) 坐标均为整张figure的归一化坐标，范围通常为0到1。
% %
% % 当前只保留一条箭头，指向s_cr临界主导分支放大框。
% % ========================================================================
% % 把原图数据坐标转换为整张图的归一化坐标，便于绘制引线箭头。
% drawnow;
% mainPosition = axesHandle.Position;
% mainToFigure = @(x,y) [ ...
%     mainPosition(1)+(x-xLower)/(xUpper-xLower)*mainPosition(3), ...
%     mainPosition(2)+(y+yLimit)/(2*yLimit)*mainPosition(4)];
% 
% % ---------- 临界分支箭尾：修改这里可以改变箭头的起始位置 ----------
% crossingSource = mainToFigure(0,0.35); % 箭尾：虚轴附近
% 
% % ---------- 临界分支箭头尖端：修改这里可以改变箭头指向 ----------
% % 目标放在左侧放大框下边缘的空白处，避免遮挡坐标和轨迹。
% crossingTarget = insetPositions(1,1:2)+[insetPositions(1,3)-0.020,-0.045];
% 
% % 如果希望临界分支箭头尖端再下移，只需减小crossingTarget第二个偏移量，
% % 例如把-0.060改成-0.080；若要向左移动，减小其第一个偏移量。
% arrowSources = crossingSource;
% arrowTargets = crossingTarget;
% arrowHandles = gobjects(1,1);
% for arrowIndex = 1
%     % 当前方向：arrowSources(arrowIndex) -> arrowTargets(arrowIndex)。
%     arrowEnd = arrowTargets(arrowIndex,:);
%     arrowVector = arrowEnd-arrowSources(arrowIndex,:);
%     % 当前箭头使用完整引线，从虚轴附近指向左侧放大框。
%     arrowStart = arrowSources(arrowIndex,:);
%     arrowHandles(arrowIndex) = annotation(figureHandle,'arrow', ...
%         [arrowStart(1),arrowEnd(1)],[arrowStart(2),arrowEnd(2)], ...
%         'Color',[0.35 0.35 0.35],'LineWidth',0.25, ...
%         'HeadStyle','vback2','HeadLength',2.5,'HeadWidth',2.5);
% end

% PNG为300dpi位图，PDF为矢量图；导出失败会直接报错而非打印成功。
% 此时主图和临界放大图均已完成，其他旧版本均在下面的块注释内。
exportgraphics(figureHandle,pngFile,'Resolution',300);
exportgraphics(figureHandle,pdfFile,'ContentType','vector');

% 已注释：上下复极点的独立双子图，当前不执行。
%{
complexFigure = figure('Color','w','Units','inches', ...
    'Position',[1 1 7.6 3.25]);
colormap(complexFigure,colorMap);
upperAxes = axes(complexFigure,'Position',[0.09 0.20 0.34 0.66]);
hold(upperAxes,'on');
plot(upperAxes,real(dominantPositive),imag(dominantPositive),'-', ...
    'Color',[0.25 0.25 0.25],'LineWidth',0.9,'HandleVisibility','off');
scatter(upperAxes,real(dominantPositive),imag(dominantPositive),24, ...
    loadingRatio,'filled','MarkerEdgeColor','none','HandleVisibility','off');
xlim(upperAxes,complexXLimits);
ylim(upperAxes,positiveImagLimits);
set(upperAxes,'FontName','Times New Roman','FontSize',10, ...
    'LineWidth',0.7,'TickDir','out','Box','on','Layer','top');
grid(upperAxes,'on');
upperAxes.GridAlpha = 0.16;
upperAxes.GridLineStyle = ':';
xlabel(upperAxes,'Real part (s^{-1})', ...
    'FontName','Times New Roman','FontSize',10);
ylabel(upperAxes,'Imaginary part (s^{-1})', ...
    'FontName','Times New Roman','FontSize',10);
title(upperAxes,'Upper-half-plane dominant pole', ...
    'FontName','Times New Roman','FontSize',10,'FontWeight','normal');

lowerAxes = axes(complexFigure,'Position',[0.52 0.20 0.34 0.66]);
hold(lowerAxes,'on');
plot(lowerAxes,real(dominantNegative),imag(dominantNegative),'-', ...
    'Color',[0.25 0.25 0.25],'LineWidth',0.9,'HandleVisibility','off');
scatter(lowerAxes,real(dominantNegative),imag(dominantNegative),24, ...
    loadingRatio,'filled','MarkerEdgeColor','none','HandleVisibility','off');
xlim(lowerAxes,complexXLimits);
ylim(lowerAxes,negativeImagLimits);
set(lowerAxes,'FontName','Times New Roman','FontSize',10, ...
    'LineWidth',0.7,'TickDir','out','Box','on','Layer','top');
grid(lowerAxes,'on');
lowerAxes.GridAlpha = 0.16;
lowerAxes.GridLineStyle = ':';
xlabel(lowerAxes,'Real part (s^{-1})', ...
    'FontName','Times New Roman','FontSize',10);
ylabel(lowerAxes,'Imaginary part (s^{-1})', ...
    'FontName','Times New Roman','FontSize',10);
title(lowerAxes,'Lower-half-plane dominant pole', ...
    'FontName','Times New Roman','FontSize',10,'FontWeight','normal');

complexColorbar = colorbar(lowerAxes);
set(complexColorbar,'Units','normalized','Position',[0.90 0.27 0.025 0.50], ...
    'FontName','Times New Roman','FontSize',10,'LineWidth',0.7);
clim(upperAxes,[0 colorLimit]);
clim(lowerAxes,[0 colorLimit]);
complexColorbar.Ticks = unique([0:0.2:colorLimit,colorLimit]);
annotation(complexFigure,'textbox',[0.87 0.79 0.13 0.07], ...
    'String','$\lambda/\lambda_4$','Interpreter','latex', ...
    'FontName','Times New Roman','FontSize',10, ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'LineStyle','none','FitBoxToText','off');

complexPngFile = fullfile(outputFolder,'dominant_root_locus_complex_zoom.png');
complexPdfFile = fullfile(outputFolder,'dominant_root_locus_complex_zoom.pdf');
try
    exportgraphics(complexFigure,complexPngFile,'Resolution',300);
    exportgraphics(complexFigure,complexPdfFile,'ContentType','vector');
catch complexExportException
    warning('The standalone complex zoom could not be exported: %s', ...
        complexExportException.message);
end
%}

% 已注释：三个分别导出的局部图，当前不执行。
%{
% Export three separate local figures in the same style as the main plot:
% the crossing point, the upper complex branch, and the lower complex branch.
localTitles = {'Virtual-axis crossing', ...
    'Upper-half-plane dominant pole', ...
    'Lower-half-plane dominant pole'};
localPngNames = {'dominant_root_locus_zoom_crossing.png', ...
    'dominant_root_locus_zoom_upper.png', ...
    'dominant_root_locus_zoom_lower.png'};
localPdfNames = {'dominant_root_locus_zoom_crossing.pdf', ...
    'dominant_root_locus_zoom_upper.pdf', ...
    'dominant_root_locus_zoom_lower.pdf'};

for localIndex = 1:3
    localFigure = figure('Color','w','Units','inches', ...
        'Position',[1 1 4.25 3.25]);
    colormap(localFigure,colorMap);
    localAxes = axes(localFigure,'Position',[0.14 0.18 0.70 0.68]);
    hold(localAxes,'on');

    switch localIndex
        case 1
            localX = zoomReal;
            localY = zeros(size(zoomReal));
            localRatio = zoomRatio;
            localXLimits = zoomXLimits;
            localYLimits = zoomYLimits;
        case 2
            localX = real(dominantPositive);
            localY = imag(dominantPositive);
            localRatio = loadingRatio;
            localXLimits = complexXLimits;
            localYLimits = positiveImagLimits;
        case 3
            localX = real(dominantNegative);
            localY = imag(dominantNegative);
            localRatio = loadingRatio;
            localXLimits = complexXLimits;
            localYLimits = negativeImagLimits;
    end

    plot(localAxes,localX,localY,'-', ...
        'Color',[0.25 0.25 0.25],'LineWidth',0.9,'HandleVisibility','off');
    if localIndex == 1
        yline(localAxes,0,'k-','LineWidth',0.60,'HandleVisibility','off');
        xline(localAxes,0,'k-','LineWidth',0.75,'HandleVisibility','off');
    end
    scatter(localAxes,localX,localY,24,localRatio, ...
        'filled','MarkerEdgeColor','none','HandleVisibility','off');
    if localIndex == 1
        plot(localAxes,0,0,'ko','MarkerSize',5,'LineWidth',0.85, ...
            'MarkerFaceColor','w','HandleVisibility','off');
    end
    xlim(localAxes,localXLimits);
    ylim(localAxes,localYLimits);
    set(localAxes,'FontName','Times New Roman','FontSize',10.5, ...
        'LineWidth',0.7,'TickDir','out','Box','on','Layer','top');
    grid(localAxes,'on');
    localAxes.GridAlpha = 0.16;
    localAxes.GridLineStyle = ':';
    xlabel(localAxes,'Real part (s^{-1})', ...
        'FontName','Times New Roman','FontSize',10);
    ylabel(localAxes,'Imaginary part (s^{-1})', ...
        'FontName','Times New Roman','FontSize',10);
    title(localAxes,localTitles{localIndex}, ...
        'FontName','Times New Roman','FontSize',10,'FontWeight','normal');

    if localIndex == 1
        text(localAxes,localXLimits(1)+0.001,localYLimits(2)-0.0002, ...
            {sprintf('$\\lambda_{\\mathrm{cr}}=%.6f$',lambdaCrossing), ...
            sprintf('$\\lambda_{\\mathrm{cr}}/\\lambda_4=%.6f$',lambdaCrossing/lambda4)}, ...
            'FontName','Times New Roman','FontSize',8, ...
            'VerticalAlignment','top','Interpreter','latex');
        text(localAxes,0.00015,0.00018,'$\\lambda_{\\mathrm{cr}}$', ...
            'FontName','Times New Roman','FontSize',8,'Interpreter','latex');
    end

    localColorbar = colorbar(localAxes);
    set(localColorbar,'Units','normalized','Position',[0.87 0.27 0.035 0.50], ...
        'FontName','Times New Roman','FontSize',10,'LineWidth',0.7);
    clim(localAxes,[0 colorLimit]);
    localColorbar.Ticks = unique([0:0.2:colorLimit,colorLimit]);
    annotation(localFigure,'textbox',[0.83 0.79 0.13 0.07], ...
        'String','$\\lambda/\\lambda_4$','Interpreter','latex', ...
        'FontName','Times New Roman','FontSize',10, ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'LineStyle','none','FitBoxToText','off');

    localPngFile = fullfile(outputFolder,localPngNames{localIndex});
    localPdfFile = fullfile(outputFolder,localPdfNames{localIndex});
    try
        exportgraphics(localFigure,localPngFile,'Resolution',300);
        exportgraphics(localFigure,localPdfFile,'ContentType','vector');
    catch localExportException
        warning('Local zoom %d could not be exported: %s',localIndex, ...
            localExportException.message);
    end
end
%}

% 已注释：四联图旧版本，当前不执行。
%{
% Combine the full locus and all local views in one four-panel figure.
combinedFigure = figure('Color','w','Units','inches', ...
    'Position',[1 1 7.6 6.2]);
colormap(combinedFigure,colorMap);
combinedLayout = tiledlayout(combinedFigure,2,2, ...
    'TileSpacing','compact','Padding','compact');

% Panel (a): full dominant-pole root locus.
fullPanel = nexttile(combinedLayout,1);
hold(fullPanel,'on');
plot(fullPanel,real(dominantPositive),imag(dominantPositive),'-', ...
    'Color',[0.25 0.25 0.25],'LineWidth',0.65,'HandleVisibility','off');
plot(fullPanel,real(dominantNegative),imag(dominantNegative),'-', ...
    'Color',[0.25 0.25 0.25],'LineWidth',0.65,'HandleVisibility','off');
plot(fullPanel,real(dominantReal),imag(dominantReal),'-', ...
    'Color',[0.25 0.25 0.25],'LineWidth',0.75,'HandleVisibility','off');
scatter(fullPanel,real(dominantPositive),imag(dominantPositive),9, ...
    loadingRatio,'filled','MarkerEdgeColor','none','HandleVisibility','off');
scatter(fullPanel,real(dominantNegative),imag(dominantNegative),9, ...
    loadingRatio,'filled','MarkerEdgeColor','none','HandleVisibility','off');
scatter(fullPanel,real(dominantReal),imag(dominantReal),13, ...
    loadingRatio,'filled','MarkerEdgeColor','none','HandleVisibility','off');
xline(fullPanel,0,'k-','LineWidth',0.65,'HandleVisibility','off');
yline(fullPanel,0,'k-','LineWidth',0.55,'HandleVisibility','off');
plot(fullPanel,0,0,'ko','MarkerSize',4,'LineWidth',0.7, ...
    'MarkerFaceColor','w','HandleVisibility','off');
xlim(fullPanel,[xLower xUpper]);
ylim(fullPanel,[-yLimit yLimit]);
set(fullPanel,'FontName','Times New Roman','FontSize',8.5, ...
    'LineWidth',0.6,'TickDir','out','Box','on','Layer','top');
grid(fullPanel,'on');
fullPanel.GridAlpha = 0.16;
fullPanel.GridLineStyle = ':';
xlabel(fullPanel,'Real part (s^{-1})','FontName','Times New Roman','FontSize',8.5);
ylabel(fullPanel,'Imaginary part (s^{-1})','FontName','Times New Roman','FontSize',8.5);
title(fullPanel,'(a) Full dominant-pole root locus', ...
    'FontName','Times New Roman','FontSize',9,'FontWeight','normal');

% Panel (b): virtual-axis crossing.
crossingPanel = nexttile(combinedLayout,2);
hold(crossingPanel,'on');
yline(crossingPanel,0,'k-','LineWidth',0.55,'HandleVisibility','off');
xline(crossingPanel,0,'k-','LineWidth',0.65,'HandleVisibility','off');
plot(crossingPanel,zoomReal,zeros(size(zoomReal)),'-', ...
    'Color',[0.25 0.25 0.25],'LineWidth',1.0,'HandleVisibility','off');
scatter(crossingPanel,zoomReal,zeros(size(zoomReal)),28,zoomRatio, ...
    'filled','MarkerEdgeColor','none','HandleVisibility','off');
plot(crossingPanel,0,0,'ko','MarkerSize',5,'LineWidth',0.85, ...
    'MarkerFaceColor','w','HandleVisibility','off');
xlim(crossingPanel,zoomXLimits);
ylim(crossingPanel,zoomYLimits);
set(crossingPanel,'FontName','Times New Roman','FontSize',8.5, ...
    'LineWidth',0.6,'TickDir','out','Box','on','Layer','top');
grid(crossingPanel,'on');
crossingPanel.GridAlpha = 0.16;
crossingPanel.GridLineStyle = ':';
xlabel(crossingPanel,'Real part (s^{-1})','FontName','Times New Roman','FontSize',8.5);
ylabel(crossingPanel,'Imaginary part (s^{-1})','FontName','Times New Roman','FontSize',8.5);
title(crossingPanel,'(b) Virtual-axis crossing', ...
    'FontName','Times New Roman','FontSize',9,'FontWeight','normal');
text(crossingPanel,zoomXLimits(1)+0.001,zoomYLimits(2)-0.0002, ...
    {sprintf('$\\lambda_{\\mathrm{cr}}=%.6f$',lambdaCrossing), ...
    sprintf('$\\lambda_{\\mathrm{cr}}/\\lambda_4=%.6f$',lambdaCrossing/lambda4)}, ...
    'FontName','Times New Roman', ...
    'FontSize',8,'VerticalAlignment','top','Interpreter','latex');
text(crossingPanel,0.00015,0.00018,'crossing', ...
    'FontName','Times New Roman','FontSize',8,'Interpreter','tex');

% Panel (c): upper complex dominant pole.
upperPanel = nexttile(combinedLayout,3);
hold(upperPanel,'on');
plot(upperPanel,real(dominantPositive),imag(dominantPositive),'-', ...
    'Color',[0.25 0.25 0.25],'LineWidth',0.9,'HandleVisibility','off');
scatter(upperPanel,real(dominantPositive),imag(dominantPositive),18, ...
    loadingRatio,'filled','MarkerEdgeColor','none','HandleVisibility','off');
xlim(upperPanel,complexXLimits);
ylim(upperPanel,positiveImagLimits);
set(upperPanel,'FontName','Times New Roman','FontSize',8.5, ...
    'LineWidth',0.6,'TickDir','out','Box','on','Layer','top');
grid(upperPanel,'on');
upperPanel.GridAlpha = 0.16;
upperPanel.GridLineStyle = ':';
xlabel(upperPanel,'Real part (s^{-1})','FontName','Times New Roman','FontSize',8.5);
ylabel(upperPanel,'Imaginary part (s^{-1})','FontName','Times New Roman','FontSize',8.5);
title(upperPanel,'(c) Upper-half-plane dominant pole', ...
    'FontName','Times New Roman','FontSize',9,'FontWeight','normal');

% Panel (d): lower complex dominant pole.
lowerPanel = nexttile(combinedLayout,4);
hold(lowerPanel,'on');
plot(lowerPanel,real(dominantNegative),imag(dominantNegative),'-', ...
    'Color',[0.25 0.25 0.25],'LineWidth',0.9,'HandleVisibility','off');
scatter(lowerPanel,real(dominantNegative),imag(dominantNegative),18, ...
    loadingRatio,'filled','MarkerEdgeColor','none','HandleVisibility','off');
xlim(lowerPanel,complexXLimits);
ylim(lowerPanel,negativeImagLimits);
set(lowerPanel,'FontName','Times New Roman','FontSize',8.5, ...
    'LineWidth',0.6,'TickDir','out','Box','on','Layer','top');
grid(lowerPanel,'on');
lowerPanel.GridAlpha = 0.16;
lowerPanel.GridLineStyle = ':';
xlabel(lowerPanel,'Real part (s^{-1})','FontName','Times New Roman','FontSize',8.5);
ylabel(lowerPanel,'Imaginary part (s^{-1})','FontName','Times New Roman','FontSize',8.5);
title(lowerPanel,'(d) Lower-half-plane dominant pole', ...
    'FontName','Times New Roman','FontSize',9,'FontWeight','normal');

clim(fullPanel,[0 colorLimit]);
clim(crossingPanel,[0 colorLimit]);
clim(upperPanel,[0 colorLimit]);
clim(lowerPanel,[0 colorLimit]);
combinedColorbar = colorbar(fullPanel);
combinedColorbar.Layout.Tile = 'east';
set(combinedColorbar,'FontName','Times New Roman','FontSize',9, ...
    'LineWidth',0.6);
combinedColorbar.Ticks = unique([0:0.2:colorLimit,colorLimit]);
title(combinedLayout,'Dominant-pole root locus and local views', ...
    'FontName','Times New Roman','FontSize',10,'FontWeight','normal');

combinedPngFile = fullfile(outputFolder,'dominant_root_locus_combined.png');
combinedPdfFile = fullfile(outputFolder,'dominant_root_locus_combined.pdf');
try
    exportgraphics(combinedFigure,combinedPngFile,'Resolution',300);
    exportgraphics(combinedFigure,combinedPdfFile,'ContentType','vector');
catch combinedExportException
    warning('The combined preview could not be exported: %s', ...
        combinedExportException.message);
end
%}

% 命令行报告临界lambda、归一化比例、限制线路和本次输出的文件。
fprintf('\nStandalone dominant root locus completed.\n');
fprintf('Output folder: %s\n',outputFolder);
fprintf('VSG buses: %d, algebraic buses: %d, reduced state order: %d\n', ...
    numel(vsgBuses),numel(algebraicBuses),2*numel(vsgBuses));
fprintf('lambda_4 = %.12f\n',lambda4);
fprintf('lambda_crossing = %.12f\n',lambdaCrossing);
fprintf('lambda_crossing/lambda_4 = %.12f\n',lambdaCrossing/lambda4);
fprintf('limiting branch = %d -> %d, boundary residual = %.3e\n', ...
    branchData(limitingEdge,1),branchData(limitingEdge,2), ...
    boundaryResidual);
fprintf('Files:\n  %s\n  %s\n  %s\n', ...
    pngFile,pdfFile, ...
    fullfile(outputFolder,'dominant_root_locus_data.csv'));

% =========================================================================
% 本文件内的局部函数：平衡点、边界、线性化、主导极点和色图计算
% =========================================================================

function [lambdaBoundary,thetaBoundary,limitingEdge,limitingSign, ...
    residualBoundary,angleErrorBoundary] = find_phase_boundary_standalone( ...
    incidence,lineWeight,varpi,referenceBus,settings)
% 搜索第一条线路相角达到pi/2的边界，并返回临界加载系数和临界平衡点。
% 输入：网络B/a、基准功率varpi、参考节点及求解容差。
% 输出：边界lambda/theta、限制线路及方向、功率残差及相角误差。

% lower为已求得平衡点的可行下界，upper为待探测加载上界。
N = size(incidence,1);
lower = 0;
upper = 1;
thetaLower = zeros(N,1);
phaseLimit = pi/2-settings.phaseAngleMargin;

% 逐步增大lambda，直到连续平衡求解失败，从而得到边界的上下包围区间。
for expansion = 1:60
    [thetaCandidate,solver] = continue_equilibrium_standalone( ...
        incidence,lineWeight,varpi,referenceBus,lower,thetaLower, ...
        upper,settings,phaseLimit);
    if solver.converged
        lower = upper;
        thetaLower = thetaCandidate;
        upper = settings.boundaryExpansionFactor*upper;
        if upper > settings.boundaryMaximumLoadingFactor
            error('The phase boundary was not bracketed.');
        end
    else
        break;
    end
end

if lower == 0
    error('No nonzero phase-cohesive equilibrium was found.');
end

% 用二分法细化可行和不可行加载区间。
for iteration = 1:settings.boundaryBisectionIterations
    middle = 0.5*(lower+upper);
    [thetaMiddle,solver] = continue_equilibrium_standalone( ...
        incidence,lineWeight,varpi,referenceBus,lower,thetaLower, ...
        middle,settings,phaseLimit);
    if solver.converged
        lower = middle;
        thetaLower = thetaMiddle;
    else
        upper = middle;
    end
end

% 相角绝对值最大的线路最可能先到pi/2，只对候选线路求增广方程。
edgeAngleLower = full(incidence.'*thetaLower);
[~,edgeOrder] = sort(abs(edgeAngleLower),'descend');
candidateEdges = edgeOrder(1:min(settings.boundaryCandidateEdges, ...
    numel(edgeOrder)));

% 用结构体保存各候选线路的解；未通过残差和相角检查的解不参与选择。
solutions = repmat(struct('valid',false,'lambda',Inf,'theta',[], ...
    'edge',NaN,'sign',NaN,'residual',Inf,'angleError',Inf), ...
    numel(candidateEdges),1);
for candidate = 1:numel(candidateEdges)
    edge = candidateEdges(candidate);
    % 相角可达到+pi/2或-pi/2，用符号把两者统一为一条边界方程。
    edgeSign = sign(edgeAngleLower(edge));
    if edgeSign == 0
        edgeSign = 1;
    end
    solutions(candidate) = solve_augmented_boundary_standalone( ...
        incidence,lineWeight,varpi,referenceBus,thetaLower,lower,upper, ...
        edge,edgeSign,settings);
end

% 选有效候选中最小lambda。若均失败，退回连续追踪给出的近边界点；
% 后续主流程仍会检查末端极点是否真的到达零，防止把失败误认成临界。
valid = find([solutions.valid]);
if isempty(valid)
    lambdaBoundary = lower;
    thetaBoundary = thetaLower;
    [~,limitingEdge] = max(abs(edgeAngleLower));
    limitingSign = sign(edgeAngleLower(limitingEdge));
    if limitingSign == 0
        limitingSign = 1;
    end
    methodBoundary = 'continuation bisection';
else
    [~,localIndex] = min([solutions(valid).lambda]);
    chosen = solutions(valid(localIndex));
    lambdaBoundary = chosen.lambda;
    thetaBoundary = chosen.theta;
    limitingEdge = chosen.edge;
    limitingSign = chosen.sign;
    methodBoundary = 'augmented edge equations';
end

% 回代方程，报告边界解的功率平衡误差和最大相角误差。
edgeAngle = full(incidence.'*thetaBoundary);
residualFull = full(incidence*(lineWeight.*sin(edgeAngle))) - ...
    lambdaBoundary*varpi;
unknown = setdiff((1:N).',referenceBus,'stable');
residualBoundary = norm(residualFull(unknown),Inf);
angleErrorBoundary = abs(max(abs(edgeAngle))-settings.boundaryAngle);

if residualBoundary > max(settings.boundaryResidualTolerance,1e-7)
    warning('Boundary residual is %.3e.',residualBoundary);
end
if angleErrorBoundary > max(settings.boundaryAngleTolerance,1e-6)
    warning(['Boundary angle error is %.3e (%s).'], ...
        angleErrorBoundary,methodBoundary);
end
end

function solution = solve_augmented_boundary_standalone( ...
    incidence,lineWeight,varpi,referenceBus,thetaInitial,lambdaLower, ...
    lambdaUpper,edge,edgeSign,settings)
% 联立节点功率平衡方程和指定线路相角等于pi/2的方程，精确求临界点。
% 未知量y=[非参考节点相角;lambda]，比普通平衡求解多一个lambda未知量。
% 对应多增加一条edgeSign*(B'*theta)_edge=pi/2约束。

N = size(incidence,1);
unknown = setdiff((1:N).',referenceBus,'stable');
y = [thetaInitial(unknown);lambdaLower];

% 默认为无效；只有通过末尾全部数值校验后才接受此线路的边界解。
solution.valid = false;
solution.lambda = Inf;
solution.theta = [];
solution.edge = edge;
solution.sign = edgeSign;
solution.residual = Inf;
solution.angleError = Inf;

for iteration = 1:100
    % 根据当前y恢复全节点相角并构造“功率平衡+线路相角约束”残差F。
    theta = zeros(N,1);
    theta(unknown) = y(1:end-1);
    lambda = y(end);
    edgeAngle = full(incidence.'*theta);
    powerResidual = full(incidence*(lineWeight.*sin(edgeAngle))) - ...
        lambda*varpi;
    F = [powerResidual(unknown); ...
        edgeSign*edgeAngle(edge)-settings.boundaryAngle];

    if norm(F,Inf) < settings.boundaryResidualTolerance
        break;
    end

    % Y=dP/dtheta；增广雅可比最后一列为dF/dlambda=-varpi，
    % 最后一行为指定支路相角对节点相角的导数。
    Y = full(incidence*spdiags(lineWeight.*cos(edgeAngle), ...
        0,numel(edgeAngle),numel(edgeAngle))*incidence.');
    J = [Y(unknown,unknown),-varpi(unknown); ...
        edgeSign*incidence(unknown,edge).',0];
    % 避免求解近奇异的增广线性系统，失败时保留solution.valid=false。
    if rcond(full(J)) < 1e-14
        return;
    end

    % 牛顿方向J*direction=-F；线搜索要求lambda>0、相角可接受、残差降低。
    direction = -(J\F);
    accepted = false;
    step = 1;
    normF = norm(F,2);
    for lineSearch = 1:settings.lineSearchMaximumIterations
        candidate = y+step*direction;
        if candidate(end) <= 0
            step = 0.5*step;
            continue;
        end
        candidateTheta = zeros(N,1);
        candidateTheta(unknown) = candidate(1:end-1);
        candidateEdgeAngle = full(incidence.'*candidateTheta);
        if max(abs(candidateEdgeAngle)) > pi/2+1e-5
            step = 0.5*step;
            continue;
        end
        candidateResidual = full(incidence* ...
            (lineWeight.*sin(candidateEdgeAngle))) - ...
            candidate(end)*varpi;
        candidateF = [candidateResidual(unknown); ...
            edgeSign*candidateEdgeAngle(edge)-settings.boundaryAngle];
        if norm(candidateF,2) < normF
            y = candidate;
            accepted = true;
            break;
        end
        step = 0.5*step;
    end
    if ~accepted
        return;
    end
end

theta = zeros(N,1);
theta(unknown) = y(1:end-1);
lambda = y(end);
edgeAngle = full(incidence.'*theta);
powerResidual = full(incidence*(lineWeight.*sin(edgeAngle))) - ...
    lambda*varpi;
residual = norm(powerResidual(unknown),Inf);
angleError = abs(edgeSign*edgeAngle(edge)-settings.boundaryAngle);

% 检查lambda仍在包围区间附近、方程残差够小且所有相角未越过容量边界。
solution.valid = isfinite(lambda) && lambda > 0 && ...
    lambda >= 0.999*lambdaLower && lambda <= 1.001*lambdaUpper && ...
    residual <= max(settings.boundaryResidualTolerance,1e-8) && ...
    angleError <= max(settings.boundaryAngleTolerance,1e-6) && ...
    max(abs(edgeAngle)) <= pi/2+1e-6;
solution.lambda = lambda;
solution.theta = theta;
solution.residual = residual;
solution.angleError = angleError;
end

function [theta,solver] = continue_equilibrium_standalone( ...
    incidence,lineWeight,varpi,referenceBus,lambdaStart,thetaStart, ...
    lambdaTarget,settings,maxAngle)
% 沿lambda方向继续追踪平衡支路，每个小步都使用牛顿法求解。
% 输入上一平衡点和目标lambda；输出目标相角以及solver收敛状态。
% 多个小步使用前一解作初值，降低直接跳到远处lambda造成求解失败的概率。

if lambdaTarget < lambdaStart
    error('The continuation target must not be smaller than its start.');
end
if lambdaTarget == lambdaStart
    theta = thetaStart;
    solver.converged = true;
    solver.reason = 'already at target';
    return;
end

% 将这段lambda区间等分，各次求解的varpi已乘上对应lambda。
numberOfSteps = max(2,settings.continuationSubsteps);
lambdaGrid = linspace(lambdaStart,lambdaTarget,numberOfSteps+1);
theta = thetaStart;
solver.converged = true;

for step = 2:numel(lambdaGrid)
    [theta,solver] = solve_equilibrium_standalone( ...
        incidence,lineWeight,lambdaGrid(step)*varpi,referenceBus, ...
        theta,settings,maxAngle);
    if ~solver.converged
        % 立即返回失败状态，由调用者决定二分缩小区间还是报错停止。
        return;
    end
end
end

function [theta,solver] = solve_equilibrium_standalone( ...
    incidence,lineWeight,varpi,referenceBus,thetaInitial,settings,maxAngle)
% 求解B*diag(a)*sin(B'*theta)=varpi，并限制线路相角不超过给定上限。
% 此函数的varpi已经含有加载系数；采用带步长减半线搜索的牛顿法。
% 固定参考节点相角为0，保留N-1个独立功率方程。

N = size(incidence,1);
unknown = setdiff((1:N).',referenceBus,'stable');
theta = thetaInitial(:);
if numel(theta) ~= N
    theta = zeros(N,1);
end
% 平移全部相角不会改变支路差值，用它统一参考相角。
theta = theta-theta(referenceBus);

solver.converged = false;
solver.reason = 'maximum iterations';
solver.residual = Inf;
solver.maximumEdgeAngle = Inf;

for iteration = 1:settings.newtonMaximumIterations
    % 计算功率方程残差和全部线路的最大相角，用两者共同判断可行收敛。
    edgeAngle = full(incidence.'*theta);
    residualFull = full(incidence*(lineWeight.*sin(edgeAngle))) - varpi;
    residual = residualFull(unknown);
    residualInf = norm(residual,Inf);
    maximumEdgeAngle = max(abs(edgeAngle));

    if residualInf < settings.newtonTolerance && ...
            maximumEdgeAngle < maxAngle
        solver.converged = true;
        solver.reason = 'converged';
        solver.residual = residualInf;
        solver.maximumEdgeAngle = maximumEdgeAngle;
        return;
    end

    Y = full(incidence*spdiags(lineWeight.*cos(edgeAngle), ...
        0,numel(edgeAngle),numel(edgeAngle))*incidence.');
    % 删除参考节点对应行列，消除功率雅可比的统一转角零方向。
    reducedY = Y(unknown,unknown);
    if rcond(reducedY) < 1e-14
        solver.reason = 'singular Newton matrix';
        break;
    end
    % 解线性方程得到牛顿步，使用左除而非显式求逆。
    direction = -(reducedY\residual);

    accepted = false;
    step = 1;
    for lineSearch = 1:settings.lineSearchMaximumIterations
        % 从完整牛顿步开始；越过相角上限或不能降低残差就把步长减半。
        candidate = theta;
        candidate(unknown) = theta(unknown)+step*direction;
        candidate(referenceBus) = 0;
        candidateEdgeAngle = full(incidence.'*candidate);
        if max(abs(candidateEdgeAngle)) >= maxAngle
            step = 0.5*step;
            continue;
        end
        candidateResidual = full(incidence* ...
            (lineWeight.*sin(candidateEdgeAngle))) - varpi;
        if norm(candidateResidual(unknown),Inf) < residualInf
            theta = candidate;
            accepted = true;
            break;
        end
        step = 0.5*step;
    end

    if ~accepted
        solver.reason = 'line search failed';
        break;
    end
end

edgeAngle = full(incidence.'*theta);
residualFull = full(incidence*(lineWeight.*sin(edgeAngle))) - varpi;
% 迭代退出后再次回代，并用最终残差更新solver，供继续追踪/边界搜索使用。
solver.residual = norm(residualFull(unknown),Inf);
solver.maximumEdgeAngle = max(abs(edgeAngle));
solver.converged = solver.residual < max(settings.newtonTolerance,1e-8) && ...
    solver.maximumEdgeAngle < maxAngle;
end

function [poles,metrics] = reduced_vsg_poles_standalone( ...
    incidence,lineWeight,varpi,theta,lambda,vsgBuses, ...
    algebraicBuses,damping,inertia)
% 对代数节点做Schur补，得到VSG节点的降阶线性化状态矩阵并计算特征值。
% 输入lambda对应的平衡相角theta；输出s=eig(A)和数值质量指标。
% 这些是平衡点附近的连续时间小信号极点，单位s^-1。

edgeAngle = full(incidence.'*theta);
powerResidual = full(incidence*(lineWeight.*sin(edgeAngle))) - ...
    lambda*varpi;
% 对线路功率a*sin(delta)求导得到a*cos(delta)，从而构造节点功率雅可比Y。
Y = full(incidence*spdiags(lineWeight.*cos(edgeAngle), ...
    0,numel(edgeAngle),numel(edgeAngle))*incidence.');

% I为动态VSG节点、L为零注入代数节点，按这两个集合分块。
% 代数小扰动满足0=YLI*dtheta_I+YLL*dtheta_L，消元后
% Yred=YII-YIL*(YLL\YLI)，不必对代数节点引入频率动态。
YII = Y(vsgBuses,vsgBuses);
YIL = Y(vsgBuses,algebraicBuses);
YLI = Y(algebraicBuses,vsgBuses);
YLL = Y(algebraicBuses,algebraicBuses);
reducedY = YII-YIL*(YLL\YLI);
reducedY = 0.5*(reducedY+reducedY.'); % 消除浮点误差引入的微小非对称。

numberOfVSG = numel(vsgBuses);
% 状态x=[dtheta_I;domega_I]：dtheta_dot=domega，
% M*domega_dot=-Yred*dtheta-D*domega，因此A=[0,I;-M^-1*Yred,-M^-1*D]。
stateMatrix = [zeros(numberOfVSG),eye(numberOfVSG); ...
    -diag(1./inertia)*reducedY, ...
    -diag(1./inertia)*diag(damping)];

poles = eig(stateMatrix); % 29个VSG对应58阶状态矩阵的全部极点。
% 记录最大线路相角、功率残差及代数分块最小特征值，便于检查数值质量。
metrics.maximumEdgeAngle = max(abs(edgeAngle));
metrics.powerFlowResidual = norm(powerResidual,Inf);
metrics.minimumYLL = min(real(eig(0.5*(YLL+YLL.'))));
end

function dominant = select_dominant_poles_standalone(poles)
% 删除统一转角对应的一个旋转零模态，再保留实部最大的主导极点。
% 采用每个lambda重新选择最右极点的定义，不是连续追踪某一个固定模态。
% 所以接近边界时主导模态可由复极点切换为实极点；两者不能伪连成同一支。

% 所有节点相角共同平移带来一个不反映稳定性的零特征值。
% 数值上删除模最小的一个极点；边界处即使有两个近零极点也只删除一个。
[~,rotationalIndex] = min(abs(poles));
poles(rotationalIndex) = [];

% 主导极点为余下极点中实部最大者；容差用于同时保留共轭对。
maximumRealPart = max(real(poles));
selectionTolerance = max(1e-8,1e-7*max(1,max(abs(real(poles)))));
selected = poles(real(poles) >= maximumRealPart-selectionTolerance);

% 没有入选的类别记为NaN；不是该类别的极点不存在，而是它当前不主导。
dominant.positive = NaN;
dominant.negative = NaN;
dominant.real = NaN;
dominant.maximumRealPart = maximumRealPart;

% 用虚部容差区分复极点和实极点，再把复极点分成上、下半平面。
complexSelected = selected(abs(imag(selected)) > 1e-8);
realSelected = selected(abs(imag(selected)) <= 1e-8);

if ~isempty(complexSelected)
    positive = complexSelected(imag(complexSelected) >= 0);
    negative = complexSelected(imag(complexSelected) < 0);
    if ~isempty(positive)
        dominant.positive = positive(1);
    end
    if ~isempty(negative)
        dominant.negative = negative(1);
    end
end
if ~isempty(realSelected)
    % 如同时存在多个入选实极点，只返回其中最右的那个。
    [~,index] = max(real(realSelected));
    dominant.real = realSelected(index);
end
end

function orderedPoles = match_poles_by_distance_standalone(previousPoles,currentPoles)
% 用最小总距离匹配相邻加载点的极点，保持同一模态的轨迹连续。
% 极点可能在复平面内交叉，因此不按实部排名，而是求全局最小距离匹配。
previousPoles = previousPoles(:);
currentPoles = currentPoles(:);
if numel(previousPoles) ~= numel(currentPoles)
    error('Adjacent pole sets must have the same size.');
end
cost = abs(previousPoles-currentPoles.').^2;
assignment = minimum_cost_assignment_standalone(cost);
orderedPoles = currentPoles(assignment);
end

function assignment = minimum_cost_assignment_standalone(cost)
% Hungarian算法：返回每一行对应的最小代价列索引。
% 该实现不依赖额外工具箱，适用于本脚本的57阶非旋转极点集合。
[numberOfRows,numberOfColumns] = size(cost);
if numberOfRows ~= numberOfColumns
    error('The pole matching cost matrix must be square.');
end
numberOfItems = numberOfRows;
u = zeros(numberOfItems+1,1);
v = zeros(numberOfItems+1,1);
p = zeros(numberOfItems+1,1);
way = zeros(numberOfItems+1,1);

for row = 1:numberOfItems
    p(1) = row;
    column0 = 1;
    minimumValues = Inf(numberOfItems+1,1);
    used = false(numberOfItems+1,1);
    while true
        used(column0) = true;
        row0 = p(column0);
        delta = Inf;
        column1 = 1;
        for column = 2:numberOfItems+1
            if ~used(column)
                reducedCost = cost(row0,column-1)-u(row0)-v(column);
                if reducedCost < minimumValues(column)
                    minimumValues(column) = reducedCost;
                    way(column) = column0;
                end
                if minimumValues(column) < delta
                    delta = minimumValues(column);
                    column1 = column;
                end
            end
        end
        for column = 1:numberOfItems+1
            if used(column)
                u(p(column)) = u(p(column))+delta;
                v(column) = v(column)-delta;
            end
        end
        column0 = column1;
        if p(column0) == 0
            break;
        end
    end
    while true
        column1 = way(column0);
        p(column0) = p(column1);
        column0 = column1;
        if column0 == 1
            break;
        end
    end
end

assignment = zeros(numberOfItems,1);
for column = 2:numberOfItems+1
    if p(column) > 0
        assignment(p(column)) = column-1;
    end
end
end

function colorMap = loading_colormap_standalone(numberOfColors)
% 构造与参考图一致的蓝-青-黄-红加载比例颜色映射。
% 五个RGB锚点经pchip平滑插值成numberOfColors行，按加载比例查表着色。

anchors = [ ...
    0.08 0.20 0.48
    0.10 0.48 0.72
    0.28 0.72 0.72
    0.90 0.70 0.22
    0.78 0.16 0.08];
position = linspace(0,1,size(anchors,1));
query = linspace(0,1,numberOfColors).';
colorMap = interp1(position,anchors,query,'pchip');
colorMap = min(max(colorMap,0),1); % 保证插值后的RGB值仍在合法区间[0,1]。
end
