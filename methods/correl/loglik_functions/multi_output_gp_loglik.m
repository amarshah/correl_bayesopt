function[ loglik ] =  multi_output_gp_loglik( X, Y, A, ls, nvar )
  
[N,D] = size(X);
 R    = size(A,1);

K = compute_kernel( X, A, ls, nvar );

C        = chol(K);
Kinv     = solve_chol(C,eye(N*R));
logdetK  = 2*sum(log(diag(C)));

loglik = -0.5*N*R*log(2*pi)-0.5*logdetK-0.5*Y'*Kinv*Y;

