% X    is NxD
% kernel params: 
% A    is RxR lower triangular, 
% ls   is RxD lengthscales
% nvar is Rx1 noise vars

function[ bigK ] = compute_kernel( X, all_difs2, A, ls, nvar )

[N,D] = size(X);
R     = size(A,1);

% compute R latent kernels K{r}
K = {};
K3 = zeros(N,N,R);
for r=1:R
    K3(:,:,r) = reshape(eye(N),N,N,1);
end

overls2 = 1./ls./ls;   % RxD
for n=1:N-1
    difs2 = reshape(all_difs2(n,n+1:N,:),N-n,D);

    r2             = difs2*overls2';  %N-n xR
    temp           = (1+sqrt(5*r2)+5/3*r2).*exp(-sqrt(5*r2));    %N-n x R
    K3(n+1:N,n,:)  = reshape(temp,N-n,1,R);  
    K3(n,n+1:N,:)  = reshape(temp,1,N-n,R);          
end

for r=1:R
    K{r} = reshape(K3(:,:,r),N,N);
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



