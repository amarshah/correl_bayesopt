% X    is NxD
% kernel params: 
% ls   is 1xD lengthscales
% nvar is noise var

function[ bigK ] = compute_kernel_small( X, A, ls, nvar )

[N,D] = size(X);
R     = size(A,1);

% compute R latent kernels K{r}
K = {};
for r=1:R
    K{r} = eye(N);
end

for n=1:N-1
    difs = bsxfun(@times,ones(N-n,1),X(n,:)) - X(n+1:N,:);     % N-n x D
    for r=1:R
        temp          = 1./ls(r,:)./ls(r,:);
        r2            = (difs.*difs)*temp';  %N-n x1
        K{r}(n+1:N,n) = (1+sqrt(5*r2)+5/3*r2).*exp(-sqrt(5*r2));  
        K{r}(n,n+1:N) = K{r}(n+1:N,n)';  
    end
end

bigK = zeros(N*R,N*R);
% NxN blocks of size RxR

for r_1=1:R
    ind1 = r_1:R:R*N;
    for r_2=r_1:R
        temp = zeros(N,N);
        if r_1==r_2
            temp = nvar(r_1)*ones(N,N);   % covs between all data for function r_1 and r_2
        end
        
        for s=1:min(r_1,r_2)
            temp = temp + A(r_1,s)*A(r_2,s)*K{s};
        end
        
        ind2 = r_2:R:R*N;
        
        bigK(ind1,ind2) = temp;
        bigK(ind2,ind1) = temp;
    end
end



