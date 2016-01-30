% Y  is NxR

function[ loglik ] =  indep_output_gp_loglik( X, Y, all_difs2, amp, ls, nvar, mean )
  
[N,D] = size(X);
 R    = length(amp);


Kall = zeros(N,N,R);
for r=1:R
    Kall(:,:,r) = reshape(nvar(r)*eye(N),N,N,1);
end

difs2 = all_difs2(1:N,1:N,:);

difs2 = difs2 + tprod(difs2,[2,1,3],triu(ones(N)),[2,1]);

overls2 = 1./ls./ls; % R x D
r2      = tprod(difs2,[1,2,-1],overls2,[3,-1]);  % N x N x R
u       = sqrt(5*r2);   
v       = (1+u+5/3*r2).*exp(-u);
Kall    = Kall + tprod(amp',[3],v,[1,2,3]); 

% for n=1:N-1
%     difs2 = reshape(all_difs2(n,n+1:N,:),N-n,D);
%     r2    = difs2*overls2';    % N-n x R
%     u     = sqrt(5*r2);        % N-n x R
%     cols  = bsxfun(@times,(1+u+5/3*r2).*exp(-u),amp);   % N-n x R
%     
%     Kall(n+1:N,n,:) = reshape(cols,N-n,1,R);
%     Kall(n,n+1:N,:) = reshape(cols,1,N-n,R);
% end

loglik = 0;
for r=1:R
    Kr        = reshape(Kall(:,:,r),N,N);
    C         = chol(Kr);
    Kinv_r    = solve_chol(C,eye(N));
    logdetKr  = 2*sum(log(diag(C)));
    Yr        = Y(:,r) - mean(r);

    loglik    = loglik - 0.5*N*log(2*pi) - 0.5*logdetKr -0.5*Yr'*Kinv_r*Yr;
end


