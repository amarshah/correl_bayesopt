% A       R  x R
% ls      RxD lengthscales
% nvar    Rx1 noise vars
% Y       NxR  

function[ A_samples, ls_samples, nvar_samples ] = sample_hypers( A, ls, nvar, X, Y, nsamples, type )
  
[N,D] = size(X);
 R    = size(A,1);

% hyper priors: 
% A         elements are N( 0 , vA  )
% logls     elements are N( mls , vls )
% lognvar   elements are N( mvar , vvar   ) 

vA   = 0.75;   mA   = 0.2;  
vls  = 0.75;   mls  = -0.5;
vvar = 0.75;   mvar = -2.75; 

logpriorA       = @(x) -0.5*sum((x-mA).*(x-mA))/vA;
logpriorlogls   = @(x) -0.5*sum((x-mls).*(x-mls))/vls;
logpriorlognvar = @(x) -0.5*sum((x-mvar).*(x-mvar))/vvar;

% first R*(R+1)/2 elt are for A,
% next  R*D elts are for logls,
% next  R elts are for lognvar,
logprior  = @(x) logpriorA(x(1:0.5*R*(R+1))) ...
                + logpriorlogls(x(0.5*R*(R+1)+1:0.5*R*(R+1)+R*D)) ...
                + logpriorlognvar(x(0.5*R*(R+1)+R*D+1:end));

if strcmp(type,'correl') 
    logpost   = @(x) correl_output_gp_loglik_vec( x, X, Y ) + logprior(x); 
elseif strcmp(type,'indep')
    logpost   = @(x)  indep_output_gp_loglik_vec( x, X, Y ) + logprior(x);     
end

initial  = A(tril(ones(R))==1);
initial  = [initial;log(ls(:));log(nvar(:))];

samples  = slicesample(initial',nsamples,'logpdf',logpost,'burnin',50,'thin',25);

  A_samples  = {};
 ls_samples  = {};
nvar_samples = {};
for ns=1:nsamples
    temp           = tril(ones(R));
    temp(temp==1)  = samples(ns,1:R*(R+1)/2);
    A_samples{ns}  = temp;
       
    ls_samples{ns} = exp(reshape(samples(ns,R*(R+1)/2+1:R*(R+1)/2+R*D),R,D));

    nvar_samples{ns} = exp(samples(ns,R*(R+1)/2+R*D+1:end));    
end

