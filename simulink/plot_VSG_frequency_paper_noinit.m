function [fig, t, f_Hz] = plot_VSG_frequency_paper_noinit(data, userOptions)
% 将 VSG30_38_Frequency 的仿真结果绘制为论文风格曲线。
% 用法：先运行 Simulink，随后运行 plot_VSG_frequency_paper_noinit
%       plot_VSG_frequency_paper_noinit(out)  % 数据在 SimulationOutput 中时
%       plot_VSG_frequency_paper_noinit(f_vsg30_38_Hz)
% 不重新仿真，不修改 SLX，不平滑、不缩放时间或频率，不减少采样点。
% 本版仅移除初始0.2lambda4文字，其余标签居中同高，并统一加粗字体。

%% 1. 常用绘图设置（只需修改这里）
% Use the Simulink model bundled in this repository.
cfg.ModelName = 'IEEE39_29VSG';
cfg.DataName = 'f_vsg30_38_Hz';
cfg.YLim = [49.90 50.40];          % 图1相同量程：改为 [49.6 51.0]
cfg.YTickStep = 0.10;              % 图1相同刻度：改为 0.20
cfg.XLim = [];                    % []：完整真实时间；也可填 [0 18.1]
cfg.XTickStep = [];               % []：自动；当前18 s图可填 3
cfg.FigureSizeCM = [16.8 8.6];     % 宽、高，单位cm
cfg.FontSize = 12;
cfg.LineWidth = 0.75;
cfg.Channels = 1:9;               % 顺序：Bus 30,31,...,38
cfg.ShowLegend = false;           % 图1无图例
cfg.ReadScheduleFromModel = true; % 只读取参数，不执行模型更新或保存
cfg.StageStart = [0 0.05 3 6 9 12 15];
cfg.StageRatio = [0.2 0.4 0.6 0.8 1.0 1.2 1.35];
cfg.OutputDir = fullfile(fileparts(mfilename('fullpath')), 'VSG_paper_figures');
cfg.FileStem = 'VSG_frequency_paper';

% 可选：用结构体覆盖上面的设置。
if nargin >= 2 && ~isempty(userOptions)
    assert(isstruct(userOptions) && isscalar(userOptions), ...
        '第二个输入必须是标量结构体。');
    keys = fieldnames(userOptions);
    for k = 1:numel(keys)
        assert(isfield(cfg, keys{k}), '未知设置：%s', keys{k});
        cfg.(keys{k}) = userOptions.(keys{k});
    end
end

%% 2. 读取原始日志：模型已设置 To Workspace，无需另接线
if nargin < 1 || isempty(data)
    [data, dataSource] = findWorkspaceLog(cfg.DataName);
else
    dataSource = '函数输入';
end
if isa(data, 'Simulink.SimulationOutput')
    assert(any(strcmp(who(data), cfg.DataName)), ...
        'SimulationOutput 中没有变量 %s。', cfg.DataName);
    data = get(data, cfg.DataName);
end
[t, f_Hz] = unpackLog(data);
assert(~isempty(t) && numel(t) > 1, '日志为空或样本不足；请先完成仿真。');
assert(all(isfinite(t)) && all(diff(t) >= 0) && t(end) > t(1), ...
    '时间必须是有效、非递减的秒数，且至少含两个不同时刻。');
assert(isreal(f_Hz), '频率数据应为实数。');
assert(all(cfg.Channels >= 1 & cfg.Channels <= size(f_Hz,2) & ...
    cfg.Channels == floor(cfg.Channels)), '所选通道超出日志范围。');
assert(numel(cfg.YLim) == 2 && diff(cfg.YLim) > 0, 'YLim 必须递增。');
assert(isscalar(cfg.YTickStep) && cfg.YTickStep > 0, 'YTickStep 必须为正数。');
assert(numel(cfg.FigureSizeCM) == 2 && all(cfg.FigureSizeCM > 0), ...
    'FigureSizeCM 必须为两个正数。');

% 已核对模型：日志输入已乘50，单位为Hz，不能再乘50。
y0 = f_Hz(1:min(100,size(f_Hz,1)), :);
y0 = y0(isfinite(y0));
if ~isempty(y0) && median(abs(y0)) < 5
    error(['输入似乎不是 Hz 频率。请使用 f_vsg30_38_Hz，', ...
        '不要直接传入 w2、w3 等标幺频率。']);
end
if any(~isfinite(f_Hz(:)))
    warning('日志包含 NaN/Inf；不填补或平滑这些点，请检查仿真是否正常结束。');
end
fprintf('数据来源：%s；%d 个样本，%d 个通道；时间 %.6g～%.6g s。\n', ...
    dataSource, numel(t), size(f_Hz,2), t(1), t(end));

%% 3. 加载阶段：默认读取当前模型，没有模型则采用上面的明确设置
[stageStart, stageRatio, scheduleSource] = readSchedule(cfg);
assert(numel(stageStart) == numel(stageRatio) && ...
    all(isfinite(stageStart)) && all(isfinite(stageRatio)) && ...
    all(diff(stageStart) >= 0), '加载时刻与倍率不匹配。');
stageStart = stageStart(:).';
stageRatio = stageRatio(:).';
fprintf('阶段标注依据：%s。\n', scheduleSource);
fprintf('请确保日志与该加载时序来自同一次设置；脚本不改动加载时刻。\n');

if isempty(cfg.XLim)
    xLim = [t(1), t(end)];
else
    xLim = cfg.XLim;
    assert(numel(xLim) == 2 && diff(xLim) > 0, 'XLim 必须递增。');
    assert(xLim(1) < t(end) && xLim(2) > t(1), 'XLim 与数据时间不重叠。');
end
xSpan = diff(xLim);
ySpan = diff(cfg.YLim);

%% 4. 白底、Times New Roman、细线和浅色网格
fig = figure('Color','w', 'Name','VSG frequency - paper style', ...
    'NumberTitle','off', 'Units','centimeters', ...
    'Position',[3 3 cfg.FigureSizeCM], 'Renderer','painters');
ax = axes('Parent',fig, 'Units','normalized', ...
    'Position',[0.12 0.17 0.855 0.79]);
hold(ax,'on');
set(ax, 'Color','w', 'Box','on', 'LineWidth',0.70, ...
    'FontName','Times New Roman', 'FontSize',cfg.FontSize, ...
    'XColor','k', 'YColor','k', 'TickDir','in', 'TickLength',[0.012 0.012], ...
    'XMinorTick','on', 'YMinorTick','on', 'Layer','top', ...
    'XGrid','on', 'YGrid','on', 'GridColor',[0.65 0.65 0.65], ...
    'GridAlpha',0.24, 'GridLineStyle','-');
xlim(ax,xLim);
ylim(ax,cfg.YLim);
if ~isempty(cfg.XTickStep)
    assert(cfg.XTickStep > 0, 'XTickStep 必须为正数。');
    xticks(ax, ceil(xLim(1)/cfg.XTickStep)*cfg.XTickStep : ...
        cfg.XTickStep : xLim(2));
end
yTicks = cfg.YLim(1):cfg.YTickStep:cfg.YLim(2);
yticks(ax,yTicks);
yticklabels(ax,arrayfun(@(v)sprintf('%g',v), yTicks, 'UniformOutput',false));

% 额定频率参考线；不参与图例。
plot(ax,xLim,[50 50], ':', 'Color',[0.48 0.48 0.48], ...
    'LineWidth',0.60, 'HandleVisibility','off');

% 阶跃线位于曲线下层；1.35lambda4用红色虚线。
for k = 2:numel(stageStart)
    if stageStart(k) > xLim(1) && stageStart(k) < xLim(2)
        color = [0.55 0.55 0.55];
        if abs(stageRatio(k)-1.35) < 1e-8
            color = [0.85 0.10 0.10];
        end
        plot(ax,[stageStart(k) stageStart(k)],cfg.YLim, '--', ...
            'Color',color, 'LineWidth',0.65, 'HandleVisibility','off');
    end
end

% 9种区分度较高的颜色；不因波形失稳而更换或隐藏某条曲线。
colors = [0.000 0.447 0.741;
          0.850 0.325 0.098;
          0.780 0.580 0.120;
          0.494 0.184 0.556;
          0.466 0.674 0.188;
          0.301 0.650 0.833;
          0.635 0.078 0.184;
          0.000 0.550 0.500;
          0.750 0.300 0.650];
h = gobjects(numel(cfg.Channels),1);
for k = 1:numel(cfg.Channels)
    c = cfg.Channels(k);
    colorIndex = mod(c-1,size(colors,1))+1;
    h(k) = plot(ax,t,f_Hz(:,c), 'Color',colors(colorIndex,:), ...
        'LineWidth',cfg.LineWidth, 'DisplayName',sprintf('Bus %d',29+c));
end
xlabel(ax,'Time (s)', 'FontName','Times New Roman', ...
    'FontSize',cfg.FontSize+1, 'FontWeight','bold');
ylabel(ax,'Frequency (Hz)', 'FontName','Times New Roman', ...
    'FontSize',cfg.FontSize+1, 'FontWeight','bold');

%% 5. 只标注加载后的六个工况；0.2lambda4初态在正文说明
% 保留全部原始数据、初始阶段及真实阶跃时刻，仅省略第一个阶段的文字。
stageEnd = [stageStart(2:end), max(t(end),xLim(2))];
for k = 2:numel(stageStart)
    left = max(stageStart(k),xLim(1));
    right = min(stageEnd(k),xLim(2));
    if right <= left
        continue;
    end
    xText = (left+right)/2;
    yText = cfg.YLim(2)-0.025*ySpan;
    if abs(stageRatio(k)-1) < 1e-8
        label = '\lambda=\lambda_4';
    else
        label = sprintf('\\lambda=%g\\lambda_4',stageRatio(k));
    end
    textColor = [0 0 0];
    if abs(stageRatio(k)-1.35) < 1e-8
        textColor = [0.85 0.10 0.10];
    end
    text(ax,xText,yText,label, 'Interpreter','tex', ...
        'FontName','Times New Roman', 'FontSize',cfg.FontSize, ...
        'FontWeight','bold', ...
        'HorizontalAlignment','center', 'VerticalAlignment','top', ...
        'Color',textColor, 'Clipping','on');
end
if cfg.ShowLegend
    legend(ax,h,'Location','southoutside','NumColumns',3, ...
        'Box','off','FontName','Times New Roman','FontSize',cfg.FontSize-1);
end

% 刻度、坐标标题、阶段标注和图例统一加粗，须在导出前设置。
set(findall(fig, '-property', 'FontWeight'), 'FontWeight', 'bold');
drawnow;

%% 6. 导出：矢量PDF、600dpi PNG、可编辑FIG和未缩放原始数据
if ~exist(cfg.OutputDir,'dir')
    [ok,msg] = mkdir(cfg.OutputDir);
    assert(ok,'无法建立输出文件夹：%s',msg);
end
base = fullfile(cfg.OutputDir,cfg.FileStem);
save([base '_data.mat'],'t','f_Hz','cfg','stageStart','stageRatio', ...
    'dataSource','scheduleSource');
savefig(fig,[base '.fig']);
drawnow;
if exist('exportgraphics','file') ~= 0
    exportgraphics(fig,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
    exportgraphics(fig,[base '.png'],'Resolution',600,'BackgroundColor','white');
else
    set(fig,'PaperUnits','centimeters','PaperPositionMode','manual', ...
        'PaperPosition',[0 0 cfg.FigureSizeCM],'PaperSize',cfg.FigureSizeCM, ...
        'InvertHardcopy','off');
    print(fig,[base '.pdf'],'-dpdf','-painters');
    print(fig,[base '.png'],'-dpng','-r600');
end
fprintf('已导出到：%s\n',cfg.OutputDir);
yShown = f_Hz(t>=xLim(1) & t<=xLim(2),cfg.Channels);
if any(yShown(:)<cfg.YLim(1) | yShown(:)>cfg.YLim(2))
    fprintf(['提示：部分曲线超出所选纵轴范围，仅在图面外被裁剪；', ...
        'FIG和MAT仍保存完整数据，未限幅。\n']);
end
end

%% 以下为数据读取和参数读取，不需要修改
function [data, source] = findWorkspaceLog(name)
if evalin('base',sprintf('exist(''%s'',''var'')',name))
    data = evalin('base',name);
    source = name;
    return;
end
candidates = {'out','simOut','simout','simulationOutput'};
for k = 1:numel(candidates)
    var = candidates{k};
    if ~evalin('base',sprintf('exist(''%s'',''var'')',var))
        continue;
    end
    obj = evalin('base',var);
    if isa(obj,'Simulink.SimulationOutput') && any(strcmp(who(obj),name))
        data = get(obj,name);
        source = [var '.' name];
        return;
    elseif isstruct(obj) && isscalar(obj) && isfield(obj,name)
        data = obj.(name);
        source = [var '.' name];
        return;
    end
end
error(['未找到 ',name,'。请先运行本次模型，不要 clear 工作区。', ...
    '\n若结果保存在其他变量中，请调用 plot_VSG_frequency_paper_noinit(你的数据)。']);
end

function [t, y] = unpackLog(data)
if isa(data,'Simulink.SimulationData.Signal')
    data = data.Values;
end
if isa(data,'timeseries')
    t = double(data.Time(:));
    y = sampleRows(data.Data,numel(t));
elseif isstruct(data) && isscalar(data) && isfield(data,'time') && isfield(data,'signals')
    t = double(data.time(:));
    assert(~isempty(t),'请使用 Timeseries 或 Structure With Time 格式保存数据。');
    cols = cell(1,numel(data.signals));
    for k = 1:numel(cols)
        cols{k} = sampleRows(data.signals(k).values,numel(t));
    end
    y = cat(2,cols{:});
else
    error('数据类型 %s 不支持；请传入 timeseries 或含time的日志结构体。',class(data));
end
end

function y = sampleRows(values, n)
y = squeeze(double(values));
if size(y,1) == n
    y = reshape(y,n,[]);
elseif ismatrix(y) && size(y,2) == n
    y = y.';
else
    error('频率数组维度与时间长度不对应。');
end
end

function [starts, ratios, source] = readSchedule(cfg)
starts = cfg.StageStart;
ratios = cfg.StageRatio;
source = '脚本的 StageStart / StageRatio';
if ~cfg.ReadScheduleFromModel
    return;
end
try
    [~,model] = fileparts(cfg.ModelName);
    if isempty(model), model = cfg.ModelName; end
    if ~bdIsLoaded(model)
        here = fileparts(mfilename('fullpath'));
        p = fullfile(here,[model '.slx']);
        if exist(p,'file')
            load_system(p);
        else
            load_system(model);
        end
    end
    sys = [model '/IEEE39_Reduced29VSG_System'];
    mode = get_param(sys,'LoadMode_common');
    if strcmp(mode,'Step sequence')
        dt = readNumber(sys,'Tstage_common');
        first = readNumber(sys,'StepTime_common');
        assert(dt>0 && first>=0 && first<dt,'无效的阶段时刻。');
        starts = [0 first (1:5)*dt];
        ratios = [0.2 0.4 0.6 0.8 1.0 1.2 1.35];
    elseif strcmp(mode,'Single step')
        starts = [0 readNumber(sys,'StepTime_common')];
        ratios = [readNumber(sys,'InitialRatio_common'),readNumber(sys,'TargetRatio_common')];
    elseif strcmp(mode,'Fixed lambda')
        starts = 0;
        ratios = readNumber(sys,'LambdaRatio_common');
    else
        error('不支持的加载方式：%s',mode);
    end
    source = [model ' 的公共加载设置'];
catch ME
    starts = cfg.StageStart;
    ratios = cfg.StageRatio;
    warning('读取模型加载参数失败，使用脚本内的时刻；请核对。原因：%s',ME.message);
end
end

function value = readNumber(block, name)
expression = get_param(block,name);
value = str2double(expression);
if ~isfinite(value)
    value = slResolve(expression,block);
end
assert(isnumeric(value) && isscalar(value) && isfinite(value), ...
    '参数 %s 不是有限数值标量。',name);
end
