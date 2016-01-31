% xstar  is 1xD
% X      is NxD
% Y      is NxR
% kernel params: 
% amp    is 1xR kernel amplitude
% ls     is RxD lengthscales
% nvar   is Rx1 noise vars
% K22invs   is NxNxR
% mean   is 1xR 

function[ m_pred, var_pred ] = compute_pred_matrices_small_noderiv( xstar, X, Y, amp, ls, K22invs, mean )

[N,D] = size(X);
 R    = size(Y,2);
 
overls2  = 1./ls./ls;    % R x D

difs     = bsxfun(@times,ones(N,1),xstar) - X;           % N x D
r2       = (difs.*difs)*overls2';              % N x R
u        = sqrt(5*r2);
K21      = bsxfun(@times,(1+u+5/3*r2).*exp(-u),amp);                    % N x R

 K22invK21  =  tprod(K22invs,[1,-1,2],K21,[-1,2]);    % N x R 

Yminusmean = Y - bsxfun(@times, ones(N,1), mean);   
    m_pred = sum(Yminusmean.*K22invK21)' + mean';  % R x 1
  var_pred = amp' - sum(K21.*K22invK21)';  % R x 1
 




     