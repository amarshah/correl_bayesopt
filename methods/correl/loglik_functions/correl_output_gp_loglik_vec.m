% Y   is NxR

function[ loglik ] = correl_output_gp_loglik_vec( x, X, Y, all_difs2 )
  
[N,D] = size(X);
 R    = size(Y,2);

% first R*(R+1)/2 elt are for A,
% next  R*D elts are for logls,
% next  R elts are for lognvar,

A       = tril(ones(R));
A(A==1) = x(1:R*(R+1)/2); 

ls      = exp(reshape(x(0.5*R*(R+1)+1:0.5*R*(R+1)+R*D),R,D));

nvar    = exp(x(0.5*R*(R+1)+R*D+1:0.5*R*(R+1)+R*D+R));

mean    = x(0.5*R*(R+1)+R*D+R+1:end);

loglik  = correl_output_gp_loglik( X, Y, all_difs2, A, ls, nvar, mean );


