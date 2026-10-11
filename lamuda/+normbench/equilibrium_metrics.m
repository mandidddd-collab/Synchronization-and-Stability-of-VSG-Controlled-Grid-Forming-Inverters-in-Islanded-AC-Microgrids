function metrics = equilibrium_metrics(data,network,lambda,theta)
%EQUILIBRIUM_METRICS Power-flow, incremental-network and reduced VSG checks.

% Compute edge-angle differences, requested injections, and electrical
% injections at the candidate equilibrium.
x = full(network.B.'*theta);
varpi = lambda*data.varpiNominal;
Pe = full(network.B*(network.a.*sin(x)));
residual = norm(Pe-varpi,Inf);
% Build the incremental weighted Laplacian evaluated at this equilibrium.
Y = full(network.B*spdiags(network.a.*cos(x),0,numel(x),numel(x))*network.B.');
eigenY = sort(real(eig(0.5*(Y+Y.'))));
if numel(eigenY) >= 2
    lambda2Y = eigenY(2);
else
    lambda2Y = NaN;
end

I = data.vsgBuses;
L = data.algebraicBuses;
YII = Y(I,I);
YIL = Y(I,L);
YLI = Y(L,I);
YLL = Y(L,L);
% The VSG and algebraic-node blocks are used for Schur-complement reduction.
minimumYLL = normbench.mineig_sym(YLL);
if minimumYLL <= 0
    % A nonpositive algebraic block cannot be eliminated reliably.
    YK = NaN(numel(I));
    reducedEigenvalues = NaN(2*numel(I),1);
    maximumRealPart = Inf;
    lambda2YK = NaN;
else
    % Eliminate algebraic angle variations to obtain the reduced VSG network.
    YK = YII-YIL*(YLL\YLI);
    YK = 0.5*(YK+YK.');
    eigYK = sort(real(eig(YK)));
    if numel(eigYK) >= 2
        lambda2YK = eigYK(2);
    else
        lambda2YK = NaN;
    end
    nI = numel(I);
    % Form the reduced first-order state matrix for the VSG dynamics.
    stateMatrix = [zeros(nI),eye(nI); ...
        -diag(1./data.M)*YK,-diag(1./data.M)*diag(data.D)];
    reducedEigenvalues = eig(stateMatrix);
    nearZeroTolerance = 1e-7;
    numberOfNearZeroModes = sum(abs(reducedEigenvalues) <= nearZeroTolerance);
    nonzero = reducedEigenvalues(abs(reducedEigenvalues) > nearZeroTolerance);
    if isempty(nonzero)
        maximumRealPart = NaN;
    else
        maximumRealPart = max(real(nonzero));
    end
end

% Return all residuals, incremental-network values, and eigenvalue data.
metrics.lambda = lambda;
metrics.theta = theta;
metrics.edgeAngle = x;
metrics.powerFlowResidual = residual;
metrics.maximumEdgeAngle = max(abs(x));
metrics.minimumEdgeCosine = min(cos(x));
metrics.Y = Y;
metrics.lambda2Y = lambda2Y;
metrics.minimumYLL = minimumYLL;
metrics.YK = YK;
metrics.lambda2YK = lambda2YK;
metrics.reducedEigenvalues = reducedEigenvalues;
metrics.maximumRealPart = maximumRealPart;
if ~exist('numberOfNearZeroModes','var')
    numberOfNearZeroModes = NaN;
end
metrics.numberOfNearZeroModes = numberOfNearZeroModes;
end
