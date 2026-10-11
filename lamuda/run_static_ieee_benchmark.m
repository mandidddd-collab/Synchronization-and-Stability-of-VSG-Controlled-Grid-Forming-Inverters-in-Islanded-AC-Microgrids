function [resultsTable,caseResults,cfg] = run_static_ieee_benchmark(cfg)
%RUN_STATIC_IEEE_BENCHMARK Run the four final tests and the independent phase boundary.
root = fileparts(mfilename('fullpath'));
addpath(root,'-begin');
if nargin < 1 || isempty(cfg), cfg = normbench.default_config(root); end
selected = normbench.prepare(cfg);
rows = repmat(struct('System',string(''),'CaseFunction',string(''), ...
    'Nodes',0,'Edges',0,'VSGBuses',0,'CPLBuses',0,'lambda_ph',NaN, ...
    'lambda_B2',NaN,'lambda_T2',NaN,'lambda_BInf',NaN,'lambda_TInf',NaN, ...
    'ratio_B2_pct',NaN,'ratio_T2_pct',NaN,'ratio_BInf_pct',NaN,'ratio_TInf_pct',NaN, ...
    'Delta2_pct',NaN,'DeltaInf_pct',NaN,'lambda_B',NaN,'lambda_P',NaN, ...
    'overall_gain_pct',NaN, ...
    'p2',NaN,'pInf',NaN,'chi2',NaN,'chiInf',NaN, ...
    'stress2',NaN,'stressInf',NaN, ...
    'Evaluation_time_s',NaN,'Boundary_residual',NaN, ...
    'Boundary_angle_error',NaN,'Limiting_from_bus',NaN, ...
    'Limiting_to_bus',NaN),numel(selected),1);
caseResults = repmat(struct('data',[],'network',[],'certificate',[], ...
    'boundary',[]),numel(selected),1);
verifiedPoints = cell(numel(selected),1);
for k = 1:numel(selected)
    data = normbench.load_convert_case(selected(k),cfg);
    network = normbench.build_network(data,cfg);
    c = normbench.compute_certificates(data,network);
    boundary = normbench.find_phase_boundary(data,network,1,cfg);
    fullResidual = norm(full(network.B*(network.a.* ...
        sin(network.B.'*boundary.theta)))-boundary.lambda*data.varpiNominal,inf);
    assert(fullResidual <= cfg.boundaryResidualTolerance, ...
        '%s: boundary residual is too large.',data.key);
    assert(boundary.angleError <= cfg.boundaryAngleTolerance, ...
        '%s: boundary angle was not refined.',data.key);
    assert(all(c.limits <= boundary.lambda+cfg.limitTolerance*max(1,boundary.lambda)), ...
        '%s: a certificate exceeds the independent boundary.',data.key);
    verifiedPoints{k} = normbench.verify_certified_points(data,network,c,cfg);
    ratios = 100*c.limits/boundary.lambda;
    edge = data.edges(boundary.limitingEdge,:);
    r.System = string(data.displayName);
    r.CaseFunction = string(data.key);
    r.Nodes = data.N;
    r.Edges = size(data.edges,1);
    r.VSGBuses = numel(data.vsgBuses);
    r.CPLBuses = numel(data.cplBuses);
    r.lambda_ph = boundary.lambda;
    r.lambda_B2 = c.lambdaB2;
    r.lambda_T2 = c.lambdaT2;
    r.lambda_BInf = c.lambdaBInf;
    r.lambda_TInf = c.lambdaTInf;
    r.ratio_B2_pct = ratios(1);
    r.ratio_T2_pct = ratios(2);
    r.ratio_BInf_pct = ratios(3);
    r.ratio_TInf_pct = ratios(4);
    r.Delta2_pct = c.gain2;
    r.DeltaInf_pct = c.gainInf;
    r.lambda_B = c.lambdaBaseline;
    r.lambda_P = c.lambdaProposed;
    r.overall_gain_pct = c.overallGain;
    r.p2 = c.p2;
    r.pInf = c.pInf;
    r.chi2 = c.chi2;
    r.chiInf = c.chiInf;
    r.stress2 = c.stress2;
    r.stressInf = c.stressInf;
    r.Evaluation_time_s = c.timeSecond;
    r.Boundary_residual = fullResidual;
    r.Boundary_angle_error = boundary.angleError;
    r.Limiting_from_bus = data.busIDs(edge(1));
    r.Limiting_to_bus = data.busIDs(edge(2));
    rows(k) = r; %#ok<AGROW>
    caseResults(k) = struct('data',data,'network',network, ...
        'certificate',c,'boundary',boundary);
    normbench.export_case_data(data,cfg);
    fprintf('%-22s  B2 %6.2f%%  T2 %6.2f%%  Binf %6.2f%%  Tinf %6.2f%%  ph %.6f\n', ...
        data.displayName,ratios,boundary.lambda);
end
resultsTable = struct2table(rows);
% Manuscript order: existing 2, existing infinity, proposed 2, proposed infinity.
resultsTable.lambda_1 = resultsTable.lambda_B2;
resultsTable.lambda_2 = resultsTable.lambda_BInf;
resultsTable.lambda_3 = resultsTable.lambda_T2;
resultsTable.lambda_4 = resultsTable.lambda_TInf;
pointVerificationTable = struct2table(vertcat(verifiedPoints{:}));
writetable(pointVerificationTable, ...
    fullfile(cfg.resultsFolder,'certified_point_verification.csv'));
write_static_outputs(resultsTable,cfg);
save(fullfile(cfg.resultsFolder,'static_benchmark_workspace.mat'), ...
    'resultsTable','caseResults','pointVerificationTable','cfg');
fprintf('\nStatic study completed: %s\n',cfg.resultsFolder);
end
