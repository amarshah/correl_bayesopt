% xstar    D  x 1
% X        N  x D
% Y        NR x 1
% A        R  x R
% ls       RxD lengthscales
% nvar     Rx1 noise vars
% mean     1xR 
% alldK12  RN x R x D

%  logEIV     scalar
% dlogEIV     Dx1
function[ logEIV ] = log_EIV_approx_noderiv( xstar, X, Y, A, ls, nvar, Z, tau, mu, K11inv, mean, difs2_star )

D = length(xstar);

[ K12, K11invK12, K11invYminusmean, ...
    Kpred, Kpredinv, logdetKpred, mpred,...
    Sigma, Sigmainv, logdetSigma, lambda,...
    Kpredinvmpred, Sigmainvlambda, SigmaKpredinv ] = compute_pred_matrices( xstar, X, A, ls, Y, tau, mu, K11inv, mean, difs2_star );


logEIV = - 0.5*( mpred' * Kpredinvmpred  + logdetKpred ) + ...
           0.5*( lambda'* Sigmainvlambda + logdetSigma ) + ...
           sum(log(Z)-0.5*mu.*mu./tau+log(2*pi*tau));
      
