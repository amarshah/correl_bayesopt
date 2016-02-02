% xstar    D  x 1
% X        N  x D
% Y        N  x R
% A        R  x R
% ls       RxD lengthscales
% nvar     Rx1 noise vars
% alldK12  RN x R x D

%  logEIV     scalar
% dlogEIV     Dx1
function[ logEIV, dlogEIV ] = log_EIV_multitask( xstar, Y, Z, tau, mu, G, Ginv, K11inv, mean, K12, dK12 )

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
C            = chol(Sigmainv);
Sigma        = solve_chol(C,eye(R));
logdetSigma  = -2*sum(log(diag(C)));           

KGpredinvmpred  = KGpredinv*mpred;     % Rx1
lambda          = Sigma*(KGpredinvmpred + mu./tau);  % Rx1
Sigmainvlambda  = KGpredinvmpred + mu./tau;  % Rx1
SigmaKGpredinv  = Sigma*KGpredinv;           % RxR

logEIV = - 0.5*( mpred' * KGpredinvmpred  + logdetKGpred ) + ...
           0.5*( lambda'* Sigmainvlambda + logdetSigma ) + ...
           sum(log(Z)-0.5*mu.*mu./tau+log(2*pi*tau));
      
dlogEIV = zeros(D,1); 

K11invdK12 = K11inv*dK12;          % N x D
dmpred     = K11invdK12'*(Y-mean);     % DxR   

dKpred     = -dK12'*K11invK12-(K12'*K11invdK12)';  % Dx1

dlogdetKGpred = R*dKpred/Kpred;   %Dx1
dlogdetSigma = -sum(sum(Sigma.*G))*dKpred;  %Dx1
%%%%%%%%%%%%%%%%%%%%%%%%%%

dmpredKGpredinv = dmpred*KGpredinv; %DxR
temp = -dKpred/Kpred/Kpred;
mpreddKGpred = tprod(temp,[1],G*mpred,[2]); %DxR
lambdadSigmainv = tprod(temp,[1],G*lambda,[2]); %DxR

dlambdaSigmainv  = tprod( dKpred,[1],1/Kpred*Ginv*1/Kpred*lambda,[2]) + ...
           (mpreddKGpred + dmpredKGpredinv);               % D x R

dlogEIV  = -0.5*(mpreddKGpred+2*dmpredKGpredinv)*mpred ...
           -0.5*dlogdetKGpred ...
           +0.5*(lambdadSigmainv+2*dlambdaSigmainv)*lambda ...
           +0.5*dlogdetSigma;   % D x 1
           
      
