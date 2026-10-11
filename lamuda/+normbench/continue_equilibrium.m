function [theta,solver] = continue_equilibrium(network,varpiDirection,referenceBus, ...
    lambdaStart,thetaStart,lambdaTarget,cfg,numberOfSteps,maxAngle)
%CONTINUE_EQUILIBRIUM Parameter continuation along lambda*varpiDirection.

% Use the configured continuation resolution unless the caller overrides it.
if nargin < 8 || isempty(numberOfSteps)
    numberOfSteps = cfg.continuationSubsteps;
end
% Use the certified phase-cohesive boundary by default.  The root-locus
% extension can pass a larger angle explicitly for exploratory points.
if nargin < 9 || isempty(maxAngle)
    maxAngle = pi/2-1e-9;
end
% The continuation path is defined for nondecreasing loading.
if lambdaTarget < lambdaStart
    error('lambdaTarget must not be smaller than lambdaStart.');
end
% No Newton solve is needed if the requested target is already reached.
if lambdaTarget == lambdaStart
    theta = thetaStart;
    solver.converged = true;
    solver.residual = 0;
    solver.iterations = 0;
    solver.maximumEdgeAngle = max(abs(full(network.B.'*thetaStart)));
    solver.reason = 'no continuation required';
    return;
end

% Split the loading interval into smaller steps so Newton receives nearby
% equilibria as initial guesses.
lambdaGrid = linspace(lambdaStart,lambdaTarget,max(2,numberOfSteps+1));
theta = thetaStart;
solver.converged = true;
for k = 2:numel(lambdaGrid)
    % Solve the next equilibrium using the previous equilibrium as the
    % initial guess and enforce the phase-cohesive angle limit.
    [theta,solver] = normbench.solve_equilibrium(network, ...
        lambdaGrid(k)*varpiDirection,referenceBus,theta,cfg,maxAngle);
    if ~solver.converged
        return;
    end
end
end
