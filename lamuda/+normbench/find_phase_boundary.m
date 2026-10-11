function boundary = find_phase_boundary(data,network,initialUpper,cfg)
%FIND_PHASE_BOUNDARY Compute the first boundary of |B^T theta|<pi/2.
%
% The function first brackets and bisects the phase-cohesive branch. It then
% solves augmented equations with candidate limiting edges fixed at pi/2.

if nargin < 3 || isempty(initialUpper) || ~isfinite(initialUpper)
    % Start with a unit loading factor when no upper guess is supplied.
    initialUpper = 1;
end
initialUpper = max(initialUpper,1e-3);
varpiDirection = data.varpiNominal;
referenceBus = data.referenceBus;
N = data.N;

boundaryTimer = tic;
% lower is known feasible; upper is intended to be infeasible.
lower = 0;
upper = initialUpper;
thetaLower = zeros(N,1);

% Expand the upper bracket until continuation fails.
for expansion = 1:60
    [thetaCandidate,solver] = normbench.continue_equilibrium(network, ...
        varpiDirection,referenceBus,lower,thetaLower,upper,cfg, ...
        cfg.continuationSubsteps);
    if solver.converged
        lower = upper;
        thetaLower = thetaCandidate;
        upper = cfg.boundaryExpansionFactor*upper;
        if upper > cfg.boundaryMaximumLoadingFactor
            error('%s boundary was not bracketed below lambda=%g.', ...
                data.displayName,cfg.boundaryMaximumLoadingFactor);
        end
    else
        break;
    end
end

if lower == 0
    % The initial upper guess was too large; shrink it until a feasible
    % phase-cohesive point is found.
    trial = initialUpper;
    for shrink = 1:40
        trial = 0.5*trial;
        [thetaCandidate,solver] = normbench.continue_equilibrium(network, ...
            varpiDirection,referenceBus,0,zeros(N,1),trial,cfg, ...
            cfg.continuationSubsteps);
        if solver.converged
            lower = trial;
            thetaLower = thetaCandidate;
            upper = 2*trial;
            break;
        end
    end
end
if lower == 0
    error('%s: no nonzero phase-cohesive equilibrium was obtained.',data.displayName);
end

% Refine the feasible/infeasible bracket by bisection.
for iteration = 1:cfg.boundaryBisectionIterations
    middle = 0.5*(lower+upper);
    [thetaMiddle,solver] = normbench.continue_equilibrium(network, ...
        varpiDirection,referenceBus,lower,thetaLower,middle,cfg,12);
    if solver.converged
        lower = middle;
        thetaLower = thetaMiddle;
    else
        upper = middle;
    end
end

edgeAngleLower = full(network.B.'*thetaLower);
% Candidate limiting branches are the edges with the largest angle at the
% last feasible continuation point.
[~,order] = sort(abs(edgeAngleLower),'descend');
numberCandidates = min(cfg.boundaryCandidateEdges,numel(order));
candidateEdges = order(1:numberCandidates);

solutions = repmat(struct('valid',false,'lambda',Inf,'theta',[], ...
    'edge',NaN,'sign',NaN,'residual',Inf,'angleError',Inf),numberCandidates,1);
for k = 1:numberCandidates
    % Refine each candidate by solving the power balance together with
    % edgeSign*angle(edge) = pi/2.
    edge = candidateEdges(k);
    edgeSign = sign(edgeAngleLower(edge));
    if edgeSign == 0
        edgeSign = 1;
    end
    solutions(k) = local_augmented_boundary(network,varpiDirection, ...
        referenceBus,thetaLower,lower,upper,edge,edgeSign,cfg);
end

valid = find([solutions.valid]);
if isempty(valid)
    % If augmented refinement fails, retain the last feasible bisection
    % point as a conservative boundary estimate.
    boundary.lambda = lower;
    boundary.theta = thetaLower;
    [~,boundary.limitingEdge] = max(abs(edgeAngleLower));
    boundary.limitingSign = sign(edgeAngleLower(boundary.limitingEdge));
    boundary.method = 'continuation/bisection';
else
    % Select the smallest valid candidate boundary.
    [~,idx] = min([solutions(valid).lambda]);
    chosen = solutions(valid(idx));
    boundary.lambda = chosen.lambda;
    boundary.theta = chosen.theta;
    boundary.limitingEdge = chosen.edge;
    boundary.limitingSign = chosen.sign;
    boundary.method = 'augmented edge-boundary equations';
end

edgeAngle = full(network.B.'*boundary.theta);
% Compute final residuals and the distance to the pi/2 angle boundary.
residualFull = full(network.B*(network.a.*sin(edgeAngle)))- ...
    boundary.lambda*varpiDirection;
unknown = setdiff((1:N).',referenceBus,'stable');
boundary.residual = norm(residualFull(unknown),Inf);
boundary.maximumEdgeAngle = max(abs(edgeAngle));
boundary.angleError = abs(boundary.maximumEdgeAngle-cfg.boundaryAngle);
boundary.lowerBracket = lower;
boundary.upperBracket = upper;
boundary.bracketWidth = upper-lower;
boundary.computationTimeSecond = toc(boundaryTimer);
boundary.candidateEdges = candidateEdges;

if boundary.residual > max(cfg.boundaryResidualTolerance,1e-7)
    warning('%s boundary residual is %.3e.',data.displayName,boundary.residual);
end
if boundary.maximumEdgeAngle > pi/2+1e-5
    error('%s refined boundary lies outside the phase-cohesive closure.', ...
        data.displayName);
end
end

function solution = local_augmented_boundary(network,varpiDirection,referenceBus, ...
    thetaInitial,lambdaLower,lambdaUpper,edge,edgeSign,cfg)
%LOCAL_AUGMENTED_BOUNDARY Solve the augmented edge-boundary equations.
N = numel(varpiDirection);
unknown = setdiff((1:N).',referenceBus,'stable');
y = [thetaInitial(unknown);lambdaLower];
solution.valid = false;
solution.lambda = Inf;
solution.theta = [];
solution.edge = edge;
solution.sign = edgeSign;
solution.residual = Inf;
solution.angleError = Inf;

for iteration = 1:80
    % Build the augmented residual: nodal power balance plus one fixed edge
    % angle at the phase-cohesive boundary.
    theta = zeros(N,1);
    theta(unknown) = y(1:end-1);
    lambda = y(end);
    x = full(network.B.'*theta);
    powerResidual = full(network.B*(network.a.*sin(x)))-lambda*varpiDirection;
    F = [powerResidual(unknown);edgeSign*x(edge)-cfg.boundaryAngle];
    if norm(F,Inf) < cfg.boundaryResidualTolerance
        break;
    end

    % Differentiate the augmented equations with respect to angles and
    % loading factor.
    Y = full(network.B*spdiags(network.a.*cos(x),0,numel(x),numel(x))*network.B.');
    J = [Y(unknown,unknown),-varpiDirection(unknown); ...
         edgeSign*full(network.B(unknown,edge)).',0];
    if rcond(J) < 1e-14
        return;
    end
    % Compute the Newton correction for the augmented system.
    direction = -(J\F);

    % Use backtracking to preserve positivity of lambda and phase cohesion.
    accepted = false;
    step = 1;
    normF = norm(F,2);
    for lineSearch = 1:35
        candidate = y+step*direction;
        if candidate(end) <= 0
            step = 0.5*step;
            continue;
        end
        thetaCandidate = zeros(N,1);
        thetaCandidate(unknown) = candidate(1:end-1);
        xCandidate = full(network.B.'*thetaCandidate);
        if max(abs(xCandidate)) > pi/2+1e-4
            step = 0.5*step;
            continue;
        end
        residualCandidate = full(network.B*(network.a.*sin(xCandidate)))- ...
            candidate(end)*varpiDirection;
        FCandidate = [residualCandidate(unknown); ...
            edgeSign*xCandidate(edge)-cfg.boundaryAngle];
        if norm(FCandidate,2) < normF
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

% Reconstruct the final full angle vector and evaluate the boundary checks.
theta = zeros(N,1);
theta(unknown) = y(1:end-1);
lambda = y(end);
x = full(network.B.'*theta);
powerResidual = full(network.B*(network.a.*sin(x)))-lambda*varpiDirection;
residual = norm(powerResidual(unknown),Inf);
angleError = abs(edgeSign*x(edge)-cfg.boundaryAngle);
valid = isfinite(lambda) && lambda > 0 && ...
    lambda >= 0.999*lambdaLower && lambda <= 1.001*lambdaUpper && ...
    residual <= max(cfg.boundaryResidualTolerance,1e-8) && ...
    angleError <= max(cfg.boundaryAngleTolerance,1e-6) && ...
    max(abs(x)) <= pi/2+1e-6;

solution.valid = valid;
solution.lambda = lambda;
solution.theta = theta;
solution.residual = residual;
solution.angleError = angleError;
end
