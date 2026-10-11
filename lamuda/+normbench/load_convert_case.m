function data = load_convert_case(spec,cfg)
%LOAD_CONVERT_CASE Convert a MATPOWER case to the paper's VSG--CPL model.
%
% Conversion rules:
%   1) use the official MATPOWER case and, when available, its solved AC PF;
%   2) retain solved voltage magnitudes but neglect R, line charging and shifts;
%   3) set a_e = V_i V_j / (x_e |tau_e|) in per unit;
%   4) form the solved net active injection Pg-Pd-Gs*V^2;
%   5) remove the AC-loss surplus proportionally from positive injections;
%   6) positive balanced injections are VSG terminals; remaining buses are
%      algebraic CPL/zero-injection nodes;
%   7) scale the balanced direction uniformly as varpi(lambda)=lambda*varpi0.

validateattributes(spec,{'struct'},{'scalar'});
% Load the local MATPOWER-format case function.
mpc = feval(spec.functionName);
solved = mpc;
powerFlowUsed = false;
powerFlowSuccess = false;
powerFlowMessage = 'case-file values';

if cfg.runACPowerFlowWhenAvailable && exist('runpf','file') == 2
    % Optionally replace stored case values with a solved AC power flow.
    try
        if exist('mpoption','file') == 2
            mpopt = mpoption('verbose',0,'out.all',0,'pf.alg','NR');
            [candidate,success] = runpf(mpc,mpopt);
        else
            [candidate,success] = runpf(mpc);
        end
        powerFlowUsed = true;
        powerFlowSuccess = logical(success);
        if success
            solved = candidate;
            powerFlowMessage = 'MATPOWER AC power flow';
        elseif cfg.requireACPowerFlowSuccess
            error('%s base AC power flow did not converge.',spec.displayName);
        else
            warning('%s AC power flow did not converge; stored case values are used.', ...
                spec.displayName);
        end
    catch powerFlowError
        if cfg.requireACPowerFlowSuccess
            rethrow(powerFlowError);
        end
        warning('%s runpf failed (%s); stored case values are used.', ...
            spec.displayName,powerFlowError.message);
        powerFlowMessage = ['runpf failed: ',powerFlowError.message];
    end
end

% MATPOWER v2 column indices used below.
BUS_I = 1; BUS_TYPE = 2; PD = 3; QD = 4; GS = 5; VM = 8; VA = 9; BASE_KV = 10;
GEN_BUS = 1; PG = 2; QG = 3; GEN_STATUS = 8; PMAX = 9;
F_BUS = 1; T_BUS = 2; BR_R = 3; BR_X = 4; BR_B = 5;
TAP = 9; SHIFT = 10; BR_STATUS = 11;

bus = solved.bus;
gen = solved.gen;
branchAll = solved.branch;
baseMVA = solved.baseMVA;
N = size(bus,1);
busIDs = bus(:,BUS_I);

% Keep only in-service branches and reject invalid reactances.
onlineBranch = branchAll(:,BR_STATUS) > 0;
branch = branchAll(onlineBranch,:);
if isempty(branch)
    error('%s has no in-service AC branches.',spec.displayName);
end
if any(branch(:,BR_X) <= cfg.minimumReactance)
    bad = find(branch(:,BR_X) <= cfg.minimumReactance,1);
    error('%s branch %d has nonpositive/near-zero reactance.',spec.displayName,bad);
end

% Convert original bus IDs to consecutive internal indices.
[foundFrom,fromBus] = ismember(branch(:,F_BUS),busIDs);
[foundTo,toBus] = ismember(branch(:,T_BUS),busIDs);
if ~all(foundFrom & foundTo)
    error('%s has branch endpoints absent from its bus table.',spec.displayName);
end

tap = branch(:,TAP);
tap(abs(tap) < 1e-12) = 1;
% Phase shifts are ignored in the lossless paper model.
shiftDegree = branch(:,SHIFT);
if any(abs(shiftDegree) > 1e-10) && cfg.ignorePhaseShifts
    warning('%s contains phase shifts; they are set to zero in the paper model.', ...
        spec.displayName);
end

voltageMagnitude = bus(:,VM);
voltageMagnitude(~isfinite(voltageMagnitude) | voltageMagnitude <= 0) = 1;
reactance = branch(:,BR_X);
% Use the lossless transfer coefficient a_e = V_i V_j/(x_e |tap_e|).
lineWeight = voltageMagnitude(fromBus).*voltageMagnitude(toBus) ./ ...
    (reactance.*abs(tap));
if any(~isfinite(lineWeight) | lineWeight <= 0)
    error('%s produced a nonpositive line-transfer coefficient.',spec.displayName);
end

% Aggregate online generation by internal bus index.
onlineGen = gen(:,GEN_STATUS) > 0;
[foundGenerator,generatorBus] = ismember(gen(onlineGen,GEN_BUS),busIDs);
if ~all(foundGenerator)
    error('%s has generator entries absent from its bus table.',spec.displayName);
end
activeGenerationPU = accumarray(generatorBus,gen(onlineGen,PG)/baseMVA, ...
    [N,1],@sum,0);
reactiveGenerationPU = accumarray(generatorBus,gen(onlineGen,QG)/baseMVA, ...
    [N,1],@sum,0);
generationCapacityPU = accumarray(generatorBus,max(gen(onlineGen,PMAX),0)/baseMVA, ...
    [N,1],@sum,0);

activeDemandPU = bus(:,PD)/baseMVA;
reactiveDemandPU = bus(:,QD)/baseMVA;
activeShuntPU = bus(:,GS).*(voltageMagnitude.^2)/baseMVA;
% Form the solved AC net active injection at each bus.
netInjectionAC = activeGenerationPU-activeDemandPU-activeShuntPU;

% Remove the positive AC-loss surplus while retaining all negative injections
% and the relative participation of positive source buses.
positive = netInjectionAC > cfg.sourceTolerance;
negative = netInjectionAC < -cfg.sourceTolerance;
sourceTotalAC = sum(netInjectionAC(positive));
demandTotal = -sum(netInjectionAC(negative));
if sourceTotalAC <= 0 || demandTotal <= 0
    error('%s does not contain both positive and negative net injections.', ...
        spec.displayName);
end
% Remove the AC-loss surplus from positive injections so the paper model is
% exactly balanced and lossless.
lossSurplus = sourceTotalAC-demandTotal;
sourceRescaling = demandTotal/sourceTotalAC;
varpiNominal = netInjectionAC;
varpiNominal(positive) = sourceRescaling*netInjectionAC(positive);
varpiNominal(abs(varpiNominal) <= cfg.sourceTolerance) = 0;

if abs(sum(varpiNominal)) > cfg.balanceTolerance*max(1,norm(varpiNominal,1))
    error('%s lossless balanced-injection construction failed.',spec.displayName);
end

% Classify buses in the converted VSG-CPL model.
vsgBuses = find(varpiNominal > cfg.sourceTolerance);
algebraicBuses = find(varpiNominal <= cfg.sourceTolerance);
cplBuses = find(varpiNominal < -cfg.sourceTolerance);
zeroInjectionBuses = find(abs(varpiNominal) <= cfg.sourceTolerance);
if isempty(vsgBuses) || isempty(cplBuses)
    error('%s conversion did not yield both VSG and CPL buses.',spec.displayName);
end

% Allocate heterogeneous damping and inertia to VSG source buses.
sourceShare = varpiNominal(vsgBuses)/sum(varpiNominal(vsgBuses));
frequencyDropRadPerSec = 2*pi*cfg.frequencyDropHzAtLambda1;
Dtotal = demandTotal/frequencyDropRadPerSec;
% Preserve the aggregate damping calibration while allowing physically
% heterogeneous VSG ratings and inertia time constants.  This avoids the
% special proportional choice D_i/M_i = constant for every terminal.
shareRange = max(sourceShare)-min(sourceShare);
if shareRange < eps
    ratingPosition = 0.5*ones(size(sourceShare));
else
    ratingPosition = (sourceShare-min(sourceShare))/shareRange;
end
dampingShape = cfg.dampingShapeRange(1) + ...
    diff(cfg.dampingShapeRange)*ratingPosition;
D = Dtotal*dampingShape/sum(dampingShape);
inertiaTimeConstant = cfg.inertiaTimeConstantRange(1) + ...
    diff(cfg.inertiaTimeConstantRange)*(1-ratingPosition);
M = D.*inertiaTimeConstant;
omegaStar = 2*pi*cfg.nominalFrequencyHz*ones(numel(vsgBuses),1);
Pstar = varpiNominal(vsgBuses);

% Use the MATPOWER slack bus as the angle reference when it is a VSG.
slackRows = find(bus(:,BUS_TYPE) == 3);
if ~isempty(slackRows) && ismember(slackRows(1),vsgBuses)
    referenceBus = slackRows(1);
else
    referenceBus = vsgBuses(1);
end

% Build audit tables containing both original and converted quantities.
busRole = repmat("Algebraic zero",N,1);
busRole(cplBuses) = "CPL";
busRole(vsgBuses) = "VSG";
busTable = table((1:N).',busIDs,bus(:,BUS_TYPE),busRole, ...
    voltageMagnitude,bus(:,VA),bus(:,BASE_KV),activeDemandPU, ...
    reactiveDemandPU,activeGenerationPU,reactiveGenerationPU, ...
    generationCapacityPU,netInjectionAC,varpiNominal, ...
    'VariableNames',{'InternalBus','OriginalBusID','MATPOWERBusType','ModelRole', ...
    'VoltageMagnitudePU','VoltageAngleDegree','BaseKV','PdPU','QdPU', ...
    'PgPU','QgPU','PmaxPU','ACNetInjectionPU','BalancedInjectionPU'});

branchTable = table((1:size(branch,1)).',branch(:,F_BUS),branch(:,T_BUS), ...
    fromBus,toBus,branch(:,BR_R),branch(:,BR_X),branch(:,BR_B),tap, ...
    shiftDegree,lineWeight, ...
    'VariableNames',{'Branch','FromBusID','ToBusID','FromInternal','ToInternal', ...
    'ResistancePU','ReactancePU','ChargingPU','TapRatio','ShiftDegree', ...
    'LosslessTransferCoefficient'});

% Pack all converted data and metadata into one structure.
data.key = spec.key;
data.caseFunction = spec.functionName;
data.displayName = spec.displayName;
data.shortName = spec.shortName;
data.sourceURL = spec.url;
data.baseMVA = baseMVA;
data.N = N;
data.busIDs = busIDs;
data.edges = [fromBus,toBus];
data.lineWeight = lineWeight;
data.voltageMagnitude = voltageMagnitude;
data.varpiNominal = varpiNominal;
data.vsgBuses = vsgBuses;
data.algebraicBuses = algebraicBuses;
data.cplBuses = cplBuses;
data.zeroInjectionBuses = zeroInjectionBuses;
data.referenceBus = referenceBus;
data.sourceShare = sourceShare;
data.D = D;
data.M = M;
data.Pstar = Pstar;
data.loadDemandNominal = -varpiNominal(algebraicBuses);
data.omegaStar = omegaStar;
data.totalDemandPU = demandTotal;
data.acLossSurplusPU = lossSurplus;
data.sourceRescaling = sourceRescaling;
data.powerFlowUsed = powerFlowUsed;
data.powerFlowSuccess = powerFlowSuccess;
data.powerFlowMessage = powerFlowMessage;
data.busTable = busTable;
data.branchTable = branchTable;
data.rawCase = mpc;
data.solvedCase = solved;
end
