function c = compute_certificates(data,network)
%COMPUTE_CERTIFICATES Compute B2, T2, Binf, and Tinf.
% B2   : ||u||_2   < g(||P||_2)
% T2   : ||u||_2   < g(sqrt((2||P||_2^2+1)/3))
% Binf : ||u||_inf < g(||P||_inf)
% Tinf : ||u||_inf < g((2||P||_inf+1)/3)
% Start a timer so the certificate calculation time can be reported.
timer = tic;
% Compute the nominal balanced injection direction in branch coordinates.
P = full(network.P);
c.uNominal = full(network.B.'*network.Ldagger*data.varpiNominal);
% Measure the loading stress under the 2-norm and infinity-norm tests.
c.stress2 = norm(c.uNominal,2);
c.stressInf = norm(c.uNominal,inf);
assert(c.stress2 > 1e-14 && c.stressInf > 1e-14,'Zero loading direction.');
% Compute the network-dependent quantities in the four inequalities.
c.p2 = max(1,norm(P,2));
c.pInf = max(1,norm(P,inf));
% Tighten the arguments of g for the two proposed conditions.
c.chi2 = sqrt((2*c.p2^2+1)/3);
c.chiInf = (2*c.pInf+1)/3;
% Evaluate the scalar certificate radius function.
c.radiusB2 = normbench.gfun(c.p2);
c.radiusT2 = normbench.gfun(c.chi2);
c.radiusBInf = normbench.gfun(c.pInf);
c.radiusTInf = normbench.gfun(c.chiInf);
% Convert each admissible radius into a critical loading factor.
c.lambdaB2 = c.radiusB2/c.stress2;
c.lambdaT2 = c.radiusT2/c.stress2;
c.lambdaBInf = c.radiusBInf/c.stressInf;
c.lambdaTInf = c.radiusTInf/c.stressInf;
c.limits = [c.lambdaB2,c.lambdaT2,c.lambdaBInf,c.lambdaTInf];
% Store baseline, proposed, and percentage-gain summaries.
c.lambdaBaseline = max(c.lambdaB2,c.lambdaBInf);
c.lambdaProposed = max(c.lambdaT2,c.lambdaTInf);
c.gain2 = 100*(c.lambdaT2/c.lambdaB2-1);
c.gainInf = 100*(c.lambdaTInf/c.lambdaBInf-1);
c.overallGain = 100*(c.lambdaProposed/c.lambdaBaseline-1);
c.timeSecond = toc(timer);
% The proposed condition must not be more restrictive than its baseline.
assert(c.lambdaB2 <= c.lambdaT2+1e-10*max(c.limits),'B2 <= T2 failed.');
assert(c.lambdaBInf <= c.lambdaTInf+1e-10*max(c.limits),'Binf <= Tinf failed.');
end
