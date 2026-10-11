function [theta,solver] = solve_equilibrium(network,varpi,referenceBus,thetaInitial,cfg,maxAngle)
%SOLVE_EQUILIBRIUM Newton solve of BA sin(B^T theta)=varpi.

if nargin < 6 || isempty(maxAngle)
    % Keep Newton iterates strictly inside the phase-cohesive boundary.
    maxAngle = pi/2-1e-9;
end
N = numel(varpi);
% Remove the rotational degree of freedom by fixing one reference bus.
unknown = setdiff((1:N).',referenceBus,'stable');
theta = thetaInitial(:);
if numel(theta) ~= N
    theta = zeros(N,1);
end
theta = theta-theta(referenceBus);

% Initialize the solver diagnostic structure.
solver.converged = false;
solver.residual = Inf;
solver.iterations = 0;
solver.maximumEdgeAngle = Inf;
solver.reason = 'maximum iterations';

for iteration = 1:cfg.newtonMaximumIterations
    % Evaluate branch angles and the nodal power-flow residual.
    edgeAngle = full(network.B.'*theta);
    residualFull = full(network.B*(network.a.*sin(edgeAngle)))-varpi;
    residual = residualFull(unknown);
    residualInf = norm(residual,Inf);
    maximumEdgeAngle = max(abs(edgeAngle));

    % Stop when both the power balance and angle constraints are satisfied.
    if residualInf < cfg.newtonTolerance && maximumEdgeAngle < maxAngle
        solver.converged = true;
        solver.residual = residualInf;
        solver.iterations = iteration-1;
        solver.maximumEdgeAngle = maximumEdgeAngle;
        solver.reason = 'converged';
        return;
    end

    % The Jacobian is the incremental weighted Laplacian.
    jacobianFull = full(network.B*spdiags(network.a.*cos(edgeAngle),0, ...
        numel(edgeAngle),numel(edgeAngle))*network.B.');
    jacobian = jacobianFull(unknown,unknown);
    if rcond(jacobian) < 1e-14
        solver.reason = 'singular Jacobian';
        break;
    end
    % Compute the Newton correction on all non-reference bus angles.
    direction = -(jacobian\residual);

    % Backtracking line search prevents a trial step from increasing the
    % residual or violating the phase-cohesive angle limit.
    accepted = false;
    step = 1;
    for lineSearch = 1:cfg.lineSearchMaximumIterations
        candidate = theta;
        candidate(unknown) = theta(unknown)+step*direction;
        candidate(referenceBus) = 0;
        candidateAngle = full(network.B.'*candidate);
        if max(abs(candidateAngle)) >= maxAngle
            step = 0.5*step;
            continue;
        end
        candidateResidual = full(network.B*(network.a.*sin(candidateAngle)))-varpi;
        if norm(candidateResidual(unknown),Inf) < residualInf
            theta = candidate;
            accepted = true;
            break;
        end
        step = 0.5*step;
    end
    % A failed line search means that this continuation point is not solved.
    if ~accepted
        solver.reason = 'line search failed';
        break;
    end
end

% Recompute final diagnostics even when the iteration limit was reached.
edgeAngle = full(network.B.'*theta);
residualFull = full(network.B*(network.a.*sin(edgeAngle)))-varpi;
solver.residual = norm(residualFull(unknown),Inf);
solver.iterations = cfg.newtonMaximumIterations;
solver.maximumEdgeAngle = max(abs(edgeAngle));
solver.converged = solver.residual < max(cfg.newtonTolerance,1e-8) && ...
    solver.maximumEdgeAngle < maxAngle;
end
