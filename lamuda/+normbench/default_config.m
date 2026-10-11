function cfg = default_config(rootFolder)
%DEFAULT_CONFIG Self-contained four-test study; no CVX or MATPOWER install.
% This function centralizes all numerical, dynamic, and plotting parameters.
if nargin == 0, rootFolder = fileparts(fileparts(mfilename('fullpath'))); end
% Store the root and default data/output folders.
cfg.rootFolder = rootFolder;
cfg.dataFolder = fullfile(rootFolder,'data_matpower81');
cfg.resultsFolder = fullfile(rootFolder,'results');
cfg.figureFolder = fullfile(cfg.resultsFolder,'figures');
cfg.caseDataFolder = fullfile(cfg.resultsFolder,'converted_case_data');
cfg.caseKeys = {'case9','case14','case_ieee30','case39','case57','case118'};
% Keep the benchmark independent of optional external power-flow software.
cfg.runACPowerFlowWhenAvailable = false;
cfg.requireACPowerFlowSuccess = false;
cfg.ignorePhaseShifts = true;
cfg.minimumReactance = 1e-8;
cfg.balanceTolerance = 1e-10;
cfg.sourceTolerance = 1e-8;
cfg.pseudoinverseToleranceFactor = 1e-12;
cfg.newtonTolerance = 1e-10;
cfg.newtonMaximumIterations = 80;
cfg.lineSearchMaximumIterations = 35;
cfg.continuationSubsteps = 20;
cfg.boundaryBisectionIterations = 50;
cfg.boundaryExpansionFactor = 1.5;
cfg.boundaryMaximumLoadingFactor = 1e4;
cfg.boundaryCandidateEdges = 8;
cfg.boundaryAngle = pi/2;
cfg.boundaryResidualTolerance = 5e-9;
cfg.boundaryAngleTolerance = 5e-7;
cfg.limitTolerance = 1e-7;
% Number of points used for static continuation curves.
cfg.staticCurvePoints = 160;
% Nominal frequency and aggregate frequency-drop calibration.
cfg.nominalFrequencyHz = 50;
cfg.frequencyDropHzAtLambda1 = 0.20;
cfg.inertiaToDampingTimeConstant = 0.50;
cfg.dampingShapeRange = [0.60 1.40];
cfg.inertiaTimeConstantRange = [0.20 1.00];
% Select the dynamic formulation and numerical solver tolerances.
cfg.dynamicSolverMode = 'reduced';
cfg.dynamicStepTime = 1.0;
cfg.dynamicEndTime = 15.0;
cfg.dynamicOutputPoints = 1501;
cfg.dynamicRelativeTolerance = 1e-8;
cfg.dynamicAbsoluteTolerance = 1e-10;
cfg.dynamicMaximumStep = 0.02;
cfg.dynamicRepresentativeCase = 'case39';
cfg.locusPoints = 101;
% Multipliers of max(lambda_T2,lambda_TInf), not absolute loading factors.
% These values are overwritten by the case-study driver when needed.
cfg.loadingMultipliers = 0.2:0.2:3.0;
cfg.loadingInitialHold = 12;
cfg.loadingHoldTime = 12;
cfg.loadingMaximumHoldTime = 36;
cfg.loadingOutputStep = 0.01;
cfg.loadingConvergenceHz = 1e-4;
cfg.loadingLossSpreadHz = 0.01;
% Figure visibility and output formats.
cfg.showFigures = true;
cfg.savePNG = true;
cfg.savePDF = true;
cfg.pngResolution = 300;
end
