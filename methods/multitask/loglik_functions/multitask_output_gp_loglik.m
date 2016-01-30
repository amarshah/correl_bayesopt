% Y  is NxR

function[ loglik ] = multitask_output_gp_loglik( X, Y, all_difs2, A, ls, nvar, mean )
  
[N,D] = size(X);
 R    = size(A,1);


difs2 = all_difs2(1:N,1:N,:);
difs2 = difs2 + tprod(difs2,[2,1,3],triu(ones(N)),[2,1]);

K       = nvar*eye(N);
overls2 = 1./ls./ls; % 1 x D
r2      = tprod(difs2,[1,2,-1],overls2',[-1]);  % N x N
u       = sqrt(5*r2);   
K       = K + (1+u+5/3*r2).*exp(-u);  %NxN

G       = A*A';  % RxR

C        = chol(K);
Kinv     = solve_chol(C,eye(N));
logdetK  = 2*sum(log(diag(C)));

D        = chol(G);
Ginv     = solve_chol(D,eye(R));
logdetG  = 2*sum(log(diag(D)));

Yvec   = (Y - mean)';   %RxN
Yvec   = Yvec(:);  %N stacked vectors of Rx1 length  (need K kron G)

KkronGinv    = kron(Kinv,Ginv);
logdetKkronG = R*logdetK + N*logdetG; 

loglik = -0.5*N*R*log(2*pi)-0.5*logdetKkronG-0.5*Yvec'*KkronGinv*Yvec;

