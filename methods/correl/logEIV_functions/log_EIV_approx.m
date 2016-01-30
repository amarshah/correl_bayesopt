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
function[ logEIV, dlogEIV ] = log_EIV_approx( xstar, X, Y, A, ls, nvar, Z, tau, mu, K11inv, mean, alldK12, difs2_star )

D = length(xstar);

[ K12, K11invK12, K11invYminusmean, ...
    Kpred, Kpredinv, logdetKpred, mpred,...
    Sigma, Sigmainv, logdetSigma, lambda,...
    Kpredinvmpred, Sigmainvlambda, SigmaKpredinv ] = compute_pred_matrices( xstar, X, A, ls, Y, tau, mu, K11inv, mean, difs2_star );


logEIV = - 0.5*( mpred' * Kpredinvmpred  + logdetKpred ) + ...
           0.5*( lambda'* Sigmainvlambda + logdetSigma ) + ...
           sum(log(Z)-0.5*mu.*mu./tau+log(2*pi*tau));
      
dlogEIV = zeros(D,1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dK12   = alldK12;

dmpred = tprod(dK12,[-1,2,1],K11invYminusmean,[-1]);      % D x R
dKpred = -tprod(dK12,[-1,2,1],K11invK12,[-1,3]) ...
         -tprod(dK12,[-1,3,1],K11invK12,[-1,2]);  % D x R x R

dKpredinv = -tprod(tprod(Kpredinv,[2,-1],dKpred,[1,-1,3]),[1,2,-1],Kpredinv,[-1,3]);  % D x R x R
dSigmainv = dKpredinv;

dlogdetKpred =  sum(sum(tprod(Kpredinv,[2,3],dKpred,[1,3,2]),3),2);  % D x 1 ?
dlogdetSigma = -sum(sum(tprod(Sigma,[2,3],dSigmainv,[1,3,2]),3),2);  % D x 1 ?

lambda        = Sigma*(Kpredinvmpred + mu./tau);

dmpredKpredinv = dmpred*Kpredinv;

dlambdaSigmainv  = tprod(tprod(SigmaKpredinv,[2,-1],dKpred,[1,-1,3]),[1,2,-1],Kpredinv'*(Kpredinvmpred+mu./tau),[-1]) + ...
           (tprod(dKpredinv,[1,2,-1],mpred,[-1]) + dmpredKpredinv);              % D x R

dlogEIV  = -0.5*(tprod(mpred,[-1],dKpredinv,[1,-1,2])+2*dmpredKpredinv)*mpred ...
           -0.5*dlogdetKpred ...
           +0.5*(tprod(lambda,[-1],dSigmainv,[1,-1,2])+2*dlambdaSigmainv)*lambda...
           +0.5*dlogdetSigma;   % D x 1
           
      


     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for d=1:D  
%     dK12   = alldK12(:,:,d);
%     
%     dmpred = dK12'*K11invY;
%     dKpred = -dK12'*K11invK12-K11invK12'*dK12;
%     
%     dKpredinv = -Kpredinv*dKpred*Kpredinv; 
%     dSigmainv = dKpredinv;
%     
%     dlogdetKpred = sum(sum(Kpredinv.*dKpred'));
%     dlogdetSigma = -sum(sum(Sigma.*dSigmainv'));
%    
%     dlambda = SigmaKpredinv*dKpred*SigmaKpredinv'*(Kpredinvmpred+mu./tau) + ...
%               Sigma*(dKpredinv*mpred+Kpredinv*dmpred);
%               
%     dlogEIV(d) = - 0.5*((mpred'*dKpredinv+ 2*dmpred'*Kpredinv)*mpred + dlogdetKpred) ...
%                  + 0.5*((lambda'*dSigmainv + 2*dlambda'*Sigmainv)*lambda + dlogdetSigma);
%              
% end
%   
%       

