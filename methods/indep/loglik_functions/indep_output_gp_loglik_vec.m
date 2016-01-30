function[ loglik ] = indep_output_gp_loglik_vec( x, X, Y, all_difs2 )
  
[N,D] = size(X);
 R    = size(Y,2);

% first R elts are for logamp,
% next  R*D elts are for logls,
% next  R elts are for lognvar,

amp     = exp(x(1:R)); 
ls      = exp(reshape(x(R+1:R+R*D),R,D));
nvar    = exp(x(R+R*D+1:2*R+R*D));
mean    = exp(x(2*R+R*D+1:end));

loglik  = indep_output_gp_loglik( X, Y, all_difs2, amp, ls, nvar, mean );


