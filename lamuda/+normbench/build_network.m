function network = build_network(data,cfg)
%BUILD_NETWORK Build incidence, Laplacian, pseudoinverse and cutset projection.

N = data.N;
d = size(data.edges,1);
% Construct the node-edge incidence matrix.  Each branch contributes +1 at
% its first bus and -1 at its second bus.
B = sparse(N,d);
for e = 1:d
    B(data.edges(e,1),e) = 1;
    B(data.edges(e,2),e) = -1;
end
% A connected network must have incidence-matrix rank N-1.
if rank(full(B),1e-10) ~= N-1
    error('%s is disconnected after conversion.',data.displayName);
end

% A is the diagonal matrix of line transfer coefficients.
A = spdiags(data.lineWeight,0,d,d);
% L is the weighted graph Laplacian.
L = B*A*B.';
% Compute the Moore-Penrose inverse of L with a scale-dependent tolerance.
pinvTol = cfg.pseudoinverseToleranceFactor*max(1,norm(full(L),2));
Ldagger = pinv(full(L),pinvTol);
% P projects branch-space vectors onto the network cut space.
P = full(B.'*Ldagger*B*A);
% U is an orthonormal basis for the same cut space.
U = orth(full(B.'));
if size(U,2) ~= N-1
    error('%s cut-space basis has an unexpected dimension.',data.displayName);
end

% Verify the projection identities before returning the network model.
projectionResidual = norm(P*P-P,'fro');
rangeResidual = norm(P*U-U,'fro');
if projectionResidual > 1e-7 || rangeResidual > 1e-7
    error('%s failed cutset-projection checks (%.3e, %.3e).', ...
        data.displayName,projectionResidual,rangeResidual);
end

% Store all matrices and numerical diagnostics in one structure.
network.B = B;
network.A = A;
network.a = data.lineWeight;
network.L = full(L);
network.Ldagger = Ldagger;
network.P = P;
network.U = U;
network.R = U*U.';
network.projectionResidual = projectionResidual;
network.rangeResidual = rangeResidual;
end
