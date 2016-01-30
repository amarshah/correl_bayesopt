% Y  is NxR

function[ loglik ] = correl_output_gp_loglik( X, Y, all_difs2, A, ls, nvar, mean )
  
[N,D] = size(X);
 R    = size(A,1);

K = compute_kernel( X, all_difs2, A, ls, nvar );

C        = cholproj(K);
Kinv     = solve_chol(C,eye(N*R));
logdetK  = 2*sum(log(diag(C)));

Yvec   = (Y - bsxfun(@times, ones(N,1), mean))';
Yvec   = Yvec(:);

loglik = -0.5*N*R*log(2*pi)-0.5*logdetK-0.5*Yvec'*Kinv*Yvec;

