% IEEE 39-bus fixed-voltage lossless network: all 29 nonzero-injection buses are VSGs.
% Source: MATPOWER 8.1 case39.
clearvars;
close all;
clc;

%% 1. Direct matrix assignments (100 MVA base; all quantities are per-unit).
lambda4 = 4.15311915166592;
sourceTolerance = 1e-8;
newtonTolerance = 1e-11;
referenceBus = 31; % Fix the reference angle at zero to remove uniform phase-shift freedom.
R = 0; % Lossless network; this is already included in the purely imaginary Y39.
V39 = [
    1.0393836; 1.0484941; 1.0307077; 1.00446; 1.0060063; 1.0082256; 0.99839728; 0.99787232;
    1.038332; 1.0178431; 1.0133858; 1.000815; 1.014923; 1.012319; 1.0161854; 1.0325203;
    1.0342365; 1.0315726; 1.0501068; 0.99101054; 1.0323192; 1.0501427; 1.0451451; 1.038001;
    1.0576827; 1.0525613; 1.0383449; 1.0503737; 1.0501149; 1.0499; 0.982; 0.9841;
    0.9972; 1.0123; 1.0494; 1.0636; 1.0275; 1.0265; 1.03;
]; % Fixed voltages, ordered by bus number.
varpi39 = [
    -0.976; 0; -3.22; -5; 0; 0; -2.338; -5.22;
    -0.065; 0; 0; -0.0853; 0; 0; -3.2; -3.29;
    0; -1.58; 0; -6.8; -2.74; 0; -2.475; -3.086;
    -2.24; -1.39; -2.81; -2.06; -2.835; 2.5; 6.68671; 6.5;
    6.32; 5.08; 6.5; 5.6; 5.4; 8.3; -1.04;
]; % Net power: positive supply, negative demand.

% Explicit IEEE 39-bus nodal admittance; unspecified entries are zero.
Y39 = complex(zeros(39));
Y39(1,[1 2 39]) = 1i*[-64.330900243309 24.330900243309003 40];
Y39(2,[1 2 3 25 30]) = 1i*[24.330900243309003 -259.42156630254556 66.2251655629139 116.27906976744185 53.90109149710281];
Y39(3,[2 3 4 18]) = 1i*[66.2251655629139 -188.36149229523767 46.948356807511736 75.18796992481204];
Y39(4,[3 4 5 14]) = 1i*[46.948356807511736 -202.59273665247298 78.125 77.51937984496124];
Y39(5,[4 5 6 8]) = 1i*[78.125 -552.026098901099 384.61538461538464 89.28571428571429];
Y39(6,[5 6 7 11 31]) = 1i*[384.61538461538464 -650.1998054324213 108.69565217391305 121.95121951219511 37.38317757009346];
Y39(7,[6 7 8]) = 1i*[108.69565217391305 -326.0869565217391 217.3913043478261];
Y39(8,[5 7 8 9]) = 1i*[89.28571428571429 217.3913043478261 -334.2252279999316 27.548209366391184];
Y39(9,[8 9 39]) = 1i*[27.548209366391184 -67.54820936639118 40];
Y39(10,[10 11 13 32]) = 1i*[-508.788215483428 232.5581395348837 232.5581395348837 46.72897196261682];
Y39(11,[6 10 11 12]) = 1i*[121.95121951219511 232.5581395348837 -377.49786479420527 22.851397362948745];
Y39(12,[11 12 13]) = 1i*[22.851397362948745 -45.43021344522614 22.851397362948745];
Y39(13,[10 12 13 14]) = 1i*[232.5581395348837 22.851397362948745 -354.5565462721092 99.00990099009901];
Y39(14,[4 13 14 15]) = 1i*[77.51937984496124 99.00990099009901 -222.61223014381602 46.08294930875576];
Y39(15,[14 15 16]) = 1i*[46.08294930875576 -152.46592803216 106.38297872340425];
Y39(16,[15 16 17 19 21 24]) = 1i*[106.38297872340425 -513.5901800650562 112.35955056179776 51.282051282051285 74.07407407407408 169.49152542372883];
Y39(17,[16 17 18 27]) = 1i*[112.35955056179776 -292.11423828208535 121.95121951219511 57.80346820809249];
Y39(18,[3 17 18]) = 1i*[75.18796992481204 121.95121951219511 -197.13918943700713];
Y39(19,[16 19 20 33]) = 1i*[51.282051282051285 -177.2843165319621 68.36204539239814 65.81545346847439];
Y39(20,[19 20 34]) = 1i*[68.36204539239814 -127.03266347451917 55.060015416804326];
Y39(21,[16 21 22]) = 1i*[74.07407407407408 -145.50264550264552 71.42857142857143];
Y39(22,[21 22 23 35]) = 1i*[71.42857142857143 -242.15568538144387 104.16666666666667 68.22445846836092];
Y39(23,[22 23 24 36]) = 1i*[104.16666666666667 -169.50280112044817 28.57142857142857 36.76470588235294];
Y39(24,[16 23 24]) = 1i*[169.49152542372883 28.57142857142857 -198.0629539951574];
Y39(25,[2 25 26 37]) = 1i*[116.27906976744185 -188.26530468393804 30.959752321981423 42.05214465937763];
Y39(26,[25 26 27 28 29]) = 1i*[30.959752321981423 -136.08400961983727 68.02721088435375 21.09704641350211 16];
Y39(27,[17 26 27]) = 1i*[57.80346820809249 68.02721088435375 -125.83067909244625];
Y39(28,[26 28 29]) = 1i*[21.09704641350211 -87.32221197641601 66.2251655629139];
Y39(29,[26 28 29 38]) = 1i*[16 66.2251655629139 -143.23890890860253 62.539086929330836];
Y39(30,[2 30]) = 1i*[53.90109149710281 -55.24861878453038];
Y39(31,[6 31]) = 1i*[37.38317757009346 -40];
Y39(32,[10 32]) = 1i*[46.72897196261682 -50];
Y39(33,[19 33]) = 1i*[65.81545346847439 -70.4225352112676];
Y39(34,[20 34]) = 1i*[55.060015416804326 -55.55555555555556];
Y39(35,[22 35]) = 1i*[68.22445846836092 -69.93006993006993];
Y39(36,[23 36]) = 1i*[36.76470588235294 -36.76470588235294];
Y39(37,[25 37]) = 1i*[42.05214465937763 -43.10344827586207];
Y39(38,[29 38]) = 1i*[62.539086929330836 -64.1025641025641];
Y39(39,[1 9 39]) = 1i*[40 40 -80];

%% 2. Determine the incidence B and line weights A from Y and V.
Y = Y39;
V = V39;
N = size(Y,1);
[toBus,fromBus,b] = find(tril(imag(Y),-1)); % Read each undirected line once.
numberOfEdges = numel(b);
B = zeros(N,numberOfEdges);
for e = 1:numberOfEdges
    B(fromBus(e),e) = 1;
    B(toBus(e),e) = -1;
end
lineWeight = V(fromBus).*V(toBus).*b; % Line weight: a_e = V_i*V_j*imag(Y_ij).
A = diag(lineWeight);

% Rescale positive net injections to remove AC-loss surplus and balance the lossless model.
varpiNominal = varpi39;
positive = varpiNominal > sourceTolerance;
negative = varpiNominal < -sourceTolerance;
demandTotalPU = -sum(varpiNominal(negative));
varpiNominal(positive) = varpiNominal(positive)* ...
    demandTotalPU/sum(varpiNominal(positive));
varpiNominal(abs(varpiNominal)<=sourceTolerance) = 0;
vsgBuses = find(abs(varpiNominal)>sourceTolerance); % Source and load buses: 29 VSGs.
algebraicBuses = find(abs(varpiNominal)<=sourceTolerance); % Ten zero-injection algebraic buses.
assert(rank(B,1e-10)==N-1 && abs(sum(varpiNominal))<1e-10);

%% 3. Fixed inertia and damping, ordered by vsgBuses.
M = diag([
    0.95304675950294; 0.9536108983754; 0.84970072736529; 0.97089134095092; 0.83044729194351;
    0.91094184021596; 0.91214350800391; 0.954253924317; 0.95126851936759; 0.9676306252391;
    0.65073679785868; 0.96582656267934; 0.96969423245907; 0.95769660762385; 0.97141212982435;
    0.96418902674761; 0.96446325793506; 0.97163962015941; 0.96394176068449; 0.96964704072291;
    0.67336300821361; 0.69727834586232; 0.71938831369835; 0.84646898587113; 0.69727834586232;
    0.79854228042407; 0.81789275505501; 0.42510619210078; 0.95509559368658;
]);
D = diag([
    1.0464330927439; 1.3801788344262; 1.6449147614291; 1.2490006953158; 1.6776349321823;
    0.91094184021596; 0.91396101960819; 1.3772042734486; 1.3905897978477; 1.1362648342662;
    1.9126252494096; 1.3087893709647; 1.2693764380121; 1.3602492758765; 1.2344253465257;
    1.1080065049794; 1.3192003343862; 1.2076542977277; 1.3229185356082; 1.2700264578416;
    1.8875694330729; 1.860029563126; 1.8334794233882; 1.65057846075; 1.860029563126;
    1.727278864437; 1.6977787091728; 2.1255309605039; 1.0559516878721;
]);
numberOfVSG = numel(vsgBuses);

%% 4. Continue equilibria, form J and select the dominant poles.
% Bus 31 has only one line, 6--31; its capacity gives the loading endpoint.
limitingEdge = find(fromBus==6 & toBus==31);
lambdaCrossing = lineWeight(limitingEdge)/abs(varpiNominal(31));
lambdaValues = unique([linspace(0,lambdaCrossing-0.02,241).'; ...
    lambdaCrossing-logspace(log10(0.02),log10(1e-6),160).';lambdaCrossing]);
numberOfLoadingPoints = numel(lambdaValues);
numberOfNonrotationalPoles = 2*numberOfVSG-1;
allNonrotationalPoles = complex(zeros(numberOfNonrotationalPoles,numberOfLoadingPoints));
dominantPositive = NaN(numberOfLoadingPoints,1);
dominantNegative = NaN(numberOfLoadingPoints,1);
dominantReal = NaN(numberOfLoadingPoints,1);
unknown = setdiff((1:N).',referenceBus);
theta = zeros(N,1);

for point = 1:numberOfLoadingPoints
    lambda = lambdaValues(point);

    % Newton power flow, initialized by the previous equilibrium.
    for iteration = 1:100
        angle = B.'*theta;
        F = B*A*sin(angle)-lambda*varpiNominal;
        if norm(F,inf)<newtonTolerance, break; end
        H = B*diag(lineWeight.*cos(angle))*B.';
        theta(unknown) = theta(unknown)-H(unknown,unknown)\F(unknown);
    end

    % Only converged, phase-feasible equilibria may be linearized.
    angle = B.'*theta;
    F = B*A*sin(angle)-lambda*varpiNominal;
    assert(norm(F,inf)<1e-8,'Power flow did not converge at lambda=%.6f.',lambda);
    assert(max(abs(angle))<=pi/2+1e-7,'Phase limit exceeded at lambda=%.6f.',lambda);

    % Eliminate zero-injection buses and form the VSG Jacobian.
    H = B*diag(lineWeight.*cos(angle))*B.';
    Hred = H(vsgBuses,vsgBuses)-H(vsgBuses,algebraicBuses)* ...
        (H(algebraicBuses,algebraicBuses)\H(algebraicBuses,vsgBuses));
    Hred = (Hred+Hred.')/2;
    J = [zeros(numberOfVSG),eye(numberOfVSG); -M\Hred,-M\D];
    poles = eig(J);
    [~,zeroMode] = min(abs(poles));
    poles(zeroMode) = []; % Remove only the uniform-phase zero mode.
    allNonrotationalPoles(:,point) = poles;

    % Dominant poles have the largest real part; keep real or conjugate branches.
    maximumRealPart = max(real(poles));
    tolerance = max(1e-8,1e-7*max(1,max(abs(real(poles)))));
    dominant = poles(real(poles)>=maximumRealPart-tolerance);
    upper = dominant(imag(dominant)>1e-8);
    lower = dominant(imag(dominant)<-1e-8);
    realPole = dominant(abs(imag(dominant))<=1e-8);
    if ~isempty(upper), dominantPositive(point) = upper(1); end
    if ~isempty(lower), dominantNegative(point) = lower(1); end
    if ~isempty(realPole), dominantReal(point) = max(real(realPole)); end
end
boundaryResidual = norm(F,inf);

%% 5. Basic root-locus plot; styling and annotations can be added in Visio.
figure('Color','w');
loadingRatio = lambdaValues/lambda4;
pointColors = repmat(loadingRatio.',numberOfNonrotationalPoles,1);
scatter(real(allNonrotationalPoles(:)),imag(allNonrotationalPoles(:)), ...
    4.5,pointColors(:),'filled'); % Color indicates loading ratio.
hold on;
dominantPoles = [dominantPositive,dominantNegative,dominantReal];
dominantColors = repmat(loadingRatio,1,3);
scatter(real(dominantPoles(:)),imag(dominantPoles(:)), ...
    12,dominantColors(:),'filled'); % Larger dots highlight the dominant modes.
colormap(turbo);
xline(0,'k--');
yline(0,'k:');
xlim([-3 0.5]); % Keep the original display range.
ylim([-30 30]);
xlabel('Real part (s^{-1})');
ylabel('Imaginary part (s^{-1})');
set(gca,'FontName','Times New Roman','FontSize',10);
grid on;
box on;

fprintf('\nlambda_4 = %.12f\n',lambda4);
fprintf('lambda_cr = %.12f; lambda_cr/lambda_4 = %.12f\n', ...
    lambdaCrossing,lambdaCrossing/lambda4);
fprintf('Limiting branch: %d -> %d; power-flow residual: %.3e\n', ...
    fromBus(limitingEdge),toBus(limitingEdge),boundaryResidual);
