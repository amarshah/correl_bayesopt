% X    is NxD
% kernel params: 
% A    is RxR lower triangular, 
% ls   is RxD lengthscales
% nvar is Rx1 noise vars

function[ bigK12, bigK22 ] = compute_kernel_partial( X, xstar, A, ls, difs2_star )

[N,D] = size(X);
R     = size(A,1);

% compute R latent kernels K{r}
K = {};
% for r=1:R
%     K{r} = eye(N);
% end
K3 = zeros(N,R);

overls2 = 1./ls./ls;   % RxD
% difs    = bsxfun(@times,ones(N,1),xstar) - X;     % N x D
% difs2   = difs.*difs;
r2      = difs2_star*overls2';  %NxR
u       = sqrt(5*r2);
K12s    = (1+u+5/3*r2).*exp(-u);    %N x R

bigK12 = zeros(N*R,R);
% NxN blocks of size RxR

for r_1=1:R
    ind1 = r_1:R:R*N;
    for r_2=1:R
        temp = zeros(N,1);
        for s=1:min(r_1,r_2)
            temp = temp + A(r_1,s)*A(r_2,s)*K12s(:,s);
        end
                
        bigK12(ind1,r_2) = temp;
    end
end

bigK22 = A*A';

