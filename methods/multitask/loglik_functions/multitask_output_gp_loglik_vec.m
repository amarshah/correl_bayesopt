% Y   is NxR

function[ loglik ] = multitask_output_gp_loglik_vec( x, X, Y, all_difs2 )
  
[N,D] = size(X);
 R    = size(Y,2);

% first R*(R+1)/2 elt are for A,
% next  D elts are for logls,
% next  1 elt  are for lognvar,

A       = tril(ones(R));
A(A==1) = x(1:R*(R+1)/2); 

ls      = exp(reshape(x(0.5*R*(R+1)+1:0.5*R*(R+1)+D),1,D));

nvar    = exp(x(end-1));

mean    = x(end);

loglik  = multitask_output_gp_loglik( X, Y, all_difs2, A, ls, nvar, mean );


