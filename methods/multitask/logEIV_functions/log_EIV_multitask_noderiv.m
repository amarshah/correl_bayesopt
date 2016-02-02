% xstar    D  x 1
% X        N  x D
% Y        N  x R
% A        R  x R
% ls       RxD lengthscales
% nvar     Rx1 noise vars
% alldK12  RN x R x D

%  logEIV     scalar
% dlogEIV     Dx1
function[ logEIV ] = log_EIV_multitask_noderiv( xstar, Y, Z, tau, mu, G, Ginv, K11inv, mean, K12 )

D = length(xstar);
R = size(G,1);

K22  = 1;

K11invK12 = K11inv*K12;          % N x 1
Kpred     = K22-K12'*K11invK12;  % scalar
KGpred    = Kpred*G;             % R x R
KGpredinv = 1/Kpred*Ginv;        % R x R
logdetKGpred = R*log(Kpred) + log(det(G));
mpred     = (Y-mean)'*K11invK12 + mean;     % Rx1   

Sigmainv     = KGpredinv + diag(1./tau);   % RxR
C            = chol(Sigmainv + 1e-5*diag(diag(Sigmainv)));
Sigma        = solve_chol(C,eye(R));
logdetSigma  = -2*sum(log(diag(C)));           

KGpredinvmpred  = KGpredinv*mpred;     % Rx1
lambda          = Sigma*(KGpredinvmpred + mu./tau);  % Rx1
Sigmainvlambda  = KGpredinvmpred + mu./tau;  % Rx1
SigmaKGpredinv  = Sigma*KGpredinv;           % RxR

logEIV = - 0.5*( mpred' * KGpredinvmpred  + logdetKGpred ) + ...
           0.5*( lambda'* Sigmainvlambda + logdetSigma ) + ...
           sum(log(Z)-0.5*mu.*mu./tau+log(2*pi*tau));
      

       