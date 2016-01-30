function[ dbigK12 ] = compute_derivs( xstar, X, A, ls )

[N,D] = size(X);
 R    = size(A,1);

dbigK12 = zeros(R*N,R,D);
%  bigK12 = zeros(R*N,R);

% compute R latent kernels K{r}
dK12 = {}; K12 = {};
for r=1:R
    dK12{r} = zeros(N,D);
%      K12{r} = zeros(N,1);
end


difs    = bsxfun(@times,ones(N,1),xstar) - X;     % NxD
difs2   = difs.*difs;
overls2 = 1./ls./ls;  % RxD
r2s     = difs2*overls2';   %NxR
us      = sqrt(r2s);
vs      = (1+us).*exp(-us); 
for r=1:R
    dr2     = 2*bsxfun(@times,difs,overls2(r,:));    %NxD
    dK12{r} = -5/6*bsxfun(@times,vs(:,r),dr2); %NxD
end

for r1=1:R
    ind1 = r1:R:R*N;
    for r2 = 1:R
        temp1 = zeros(N,D);
%         temp2 = zeros(N,1);
%         if r1==r2
%             temp2 = nvar(r1)*ones(N,1);    % covs between xstar and all observed for function r1 and r2 
%         end
        
        for s=1:min(r1,r2)
            aa    = A(r1,s)*A(r2,s);
            temp1 = temp1 + aa*dK12{s};   % NxD
%             temp2 = temp2 + aa* K12{s};   % Nx1
        end                
        
        dbigK12(ind1,r2,:) = reshape(temp1,N,1,D);        
%          bigK12(ind1,r2) = temp2;
    end
end





    