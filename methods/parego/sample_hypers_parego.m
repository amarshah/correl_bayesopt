% amp     scalar
% ls      1xD lengthscales
% nvar    scalar
% Y       Nx1  

function[ amp_samples, ls_samples, nvar_samples, K22inv_samples, mean_samples ] = sample_hypers_parego( X, all_difs2, Y, nsamples, lower, upper) 
  
[N,D] = size(X);
 
% hyper priors: 
% logamp    elements are N( mamp , vamp  )
% logls     elements are N( mls , vls )
% lognvar   elements are N( mvar , vvar   ) 
% mean      elements are N( 0, meanvar )

vamp = 0.75;   mamp = 0.2;  
vls  = 0.3;   mls  = -0.6;
vvar = 0.5;   mvar = -1.75; 
meanmean = sum(upper+lower);  meanvar = sum((upper-lower).^2/20);

logpriorlogamp  = @(x) -0.5*sum((x-mamp).*(x-mamp))/vamp  - sum(x); % '-x' is for the log change of var
logpriorlogls   = @(x) -0.5*sum((x-mls).*(x-mls))/vls  - sum(x);
logpriorlognvar = @(x) -0.5*sum((x-mvar).*(x-mvar))/vvar  - sum(x);
logpriormean    = @(x) -0.5*sum((x-meanmean).*(x-meanmean))/meanvar;

% first 1 elt is for logamp,
% next  D elts are for logls,
% next  1 elt is for lognvar,
% next  1 elt is for mean,
logprior  = @(x) logpriorlogamp(x(1)) ...
                + logpriorlogls(x(2:D+1)) ...
                + logpriorlognvar(x(D+2)) ...
                + logpriormean(x(D+3));

logpost   = @(x) indep_output_gp_loglik_vec( x, X, Y, all_difs2 ) + logprior(x); 

initial  = [-0.25;zeros(D,1);-2.5;0];

samples  = slicesample(initial',nsamples,'logpdf',logpost,'burnin',25,'thin',10);

amp_samples  = {};
 ls_samples  = {};
nvar_samples = {};
mean_samples = {};
K22inv_samples = {};

for ns=1:nsamples
    amp_samples{ns} = exp(samples(ns,1));
       
    ls_samples{ns}  = exp(samples(ns,2:D+1));

    nvar_samples{ns} = max(exp(samples(ns,D+2)), 1e-6);    
    
    mean_samples{ns} = samples(ns,D+3);    
    
    amp  =  amp_samples{ns};
    ls   =   ls_samples{ns};
    nvar = nvar_samples{ns};

    K = (amp+nvar)*eye(N);

    overls2 = 1./ls./ls;   % 1 x D
    for n=1:N-1
        difs2         = reshape(all_difs2(n,n+1:N,:),N-n,D);
        r2            = difs2*overls2';  %N-n x1
        u             = sqrt(5*r2);
        K(n+1:N,n)    = amp*(1+u+5/3*r2).*exp(-u);  
        K(n,n+1:N)    = K(n+1:N,n)'; 
    end

    C       = chol(K);
    K22inv  = solve_chol(C,eye(N));
        
    K22inv_samples{ns} = K22inv;
    
end

        
end




