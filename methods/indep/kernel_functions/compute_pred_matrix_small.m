% xstar  is 1xD
% X      is NxD
% Y      is Nx1
% kernel params: 
% amp    is kernel amplitude
% ls     is 1xD lengthscales
% nvar   is noise var

function[ m_pred, var_pred, dm_pred, dvar_pred ] = compute_pred_matrix_small( xstar, X, Y, amp, ls, K22inv )

[N,D] = size(X);

K11 = amp;

overls2  = 1./ls./ls;    % 1 x D

difs     = bsxfun(@times,ones(N,1),xstar) - X;           % N x D
r2       = (difs.*difs)*overls2';                        % N x 1
u        = sqrt(5*r2);
K21      = amp*(1+u+5/3*r2).*exp(-u);                    % N x 1

dr2      = 2*bsxfun(@times,difs,overls2);                % N x D
r2_rep   = bsxfun(@times,r2,ones(1,D));                  % N x D
u        = sqrt(5*r2_rep);
dK21     = -5/6*dr2.*(1+u).*exp(-u);                     % N x D

% for n=1:N
%     difs = bsxfun(@times,ones(N+1-n,1),Xfull(n,:)) - Xfull(n+1:N+1,:);     % N+1-n x D
%     r2           = (difs.*difs)*temp';  %N+1-n x1
%     K(n+1:N+1,n) = amp*(1+sqrt(5*r2)+5/3*r2).*exp(-sqrt(5*r2));  
%     K(n,n+1:N+1) = K(n+1:N+1,n)';  
%     if n==1
%         dr2     = 2*difs.*bsxfun(@times,ones(N,1),temp);       %NxD
%         r2_rep  = bsxfun(@times,r2,ones(1,D));
%         dK12    = -5/6*dr2.*(1+sqrt(5*r2_rep)).*exp(-sqrt(5*r2_rep)); %NxD
%     end
% end
% 
% 
% C   = chol(K22);
% K22inv = solve_chol(C,eye(N));

 K22invK21 =  K22inv*K21;     % Nx1
K22invdK21 =  K22inv*dK21;    % NxD

   m_pred  = Y'*K22invK21;      % scalar
 var_pred  = K11 - K21'*K22invK21;  %scalar

  dm_pred  = K22invdK21'*Y;    % Dx1
dvar_pred  = - (K21'*K22invdK21)' - K22invdK21'*K21;  % Dx1




     