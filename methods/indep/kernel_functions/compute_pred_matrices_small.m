% xstar  is 1xD
% X      is NxD
% Y      is NxR
% kernel params: 
% amp    is 1xR kernel amplitude
% ls     is RxD lengthscales
% nvar   is Rx1 noise vars
% K22invs   is NxNxR
% mean   is 1xR 

function[ m_pred, var_pred, dm_pred, dvar_pred ] = compute_pred_matrices_small( xstar, X, Y, amp, ls, K22invs, mean )

[N,D] = size(X);
 R    = size(Y,2);
 
overls2  = 1./ls./ls;    % R x D

difs     = bsxfun(@times,ones(N,1),xstar) - X;           % N x D
r2       = (difs.*difs)*overls2';              % N x R
u        = sqrt(5*r2);
K21      = bsxfun(@times,(1+u+5/3*r2).*exp(-u),amp);                    % N x R

dr2      = 2*tprod(difs,[1,2],overls2,[3,2]);            % N x D x R
v        = (1+u).*exp(-u);                                % N x R
dK21     = -5/6*tprod(dr2,[1,2,3],v,[1,3]);              % N x D x R


 K22invK21  =  tprod(K22invs,[1,-1,2],K21,[-1,2]);    % N x R 
 K22invdK21 =  tprod(K22invs,[1,-1,3],dK21,[-1,2,3]);    % N x D x R  

Yminusmean = Y - bsxfun(@times, ones(N,1), mean);   
    m_pred = sum(Yminusmean.*K22invK21)' + mean';  % R x 1
  var_pred = amp' - sum(K21.*K22invK21)';  % R x 1
 
  dm_pred  = tprod(K22invdK21,[-1,1,2],Yminusmean,[-1,2]);   % D x R 
dvar_pred  = - tprod(K22invdK21,[-1,1,2],K21,[-1,2]) ...
             - tprod(tprod(K22invs,[-1,1,3],dK21,[-1,2,3]),[-1,1,2],K21,[-1,2]);    % D x R

  var_pred = max(var_pred, 1e-8);
%(K21'*K22invdK21)' - K22invdK21'*K21;  % Dx1




     