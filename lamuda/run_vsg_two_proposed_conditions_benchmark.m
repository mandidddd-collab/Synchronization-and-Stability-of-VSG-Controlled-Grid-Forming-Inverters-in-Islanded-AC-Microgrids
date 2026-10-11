function [resultsTable,caseResults,cfg] = run_vsg_two_proposed_conditions_benchmark(cfg)
%RUN_VSG_TWO_PROPOSED_CONDITIONS_BENCHMARK
% Reproducible six-system benchmark for the final four conditions:
%
%   B2   : ||u||_2   < g(||P||_2)                         [literature]
%   T2   : ||u||_2   < g(sqrt((2||P||_2^2+1)/3))          [proposed]
%   Binf : ||u||_inf < g(||P||_inf)                       [literature]
%   Tinf : ||u||_inf < g((2||P||_inf+1)/3)                 [proposed]
%
% The bundled MATPOWER case files are used as data files.  A separate
% MATPOWER installation, CVX, and SDPT3 are not required.  The direct
% phase-cohesive boundary is found
% independently by continuation and boundary refinement.
%
% Outputs:
%   results/static_benchmark_results.csv
%   results/static_benchmark_workspace.mat

root = fileparts(mfilename('fullpath'));
addpath(root,'-begin');
if nargin < 1 || isempty(cfg)
    cfg = normbench.default_config(root);
end
[resultsTable,caseResults,cfg] = run_static_ieee_benchmark(cfg);
end
