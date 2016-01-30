% amp     R  x 1
% ls      RxD lengthscales
% nvar    Rx1 noise vars
% Y       NxR  

function[ A_samples, ls_samples, nvar_samples, Kinv_samples, mean_samples, obj_mean_pred ] = sample_hypers_multitask( A, ls, nvar, mean, X, all_difs2, Y, nsamples, lower, upper )
  
[N,D] = size(X);
 R    = size(A,1);
 
% hyper priors:
% A         elements are N( 0 , vA  )
% logls     elements are N( mls , vls )
% lognvar   elements are N( mvar , vvar   ) 
% mean      elements are N( meanmean, meanvar )

vA   = 0.5;   mA   = 0; 
vls  = 0.3;   mls  = -0.6;
mvar = log(sum((upper-lower).^2/40)/R/50); vvar = abs(mvar)/2; 
meanmean = sum((upper+lower)/2)/R;  meanvar = sum((upper-lower).^2/40)/R;

logpriorA       = @(x) -0.5*sum((x-mA).*(x-mA))/vA;
logpriorlogls   = @(x) -0.5*sum((x-mls).*(x-mls))/vls  - sum(x); % '-x' is for the log change of var
logpriorlognvar = @(x) -0.5*sum((x-mvar).*(x-mvar))/vvar - sum(x);
logpriormean    = @(x) -0.5*sum((x-meanmean).*(x-meanmean))/meanvar;

% first R*(R+1)/2 elt are for A,
% next  D elts are for logls,
% next  1 elt  are for lognvar,
% next  1 elt  are for meanvar
logprior  = @(x) logpriorA(x(1:0.5*R*(R+1))) ...
                + logpriorlogls(x(0.5*R*(R+1)+1:0.5*R*(R+1)+D)) ...
                + logpriorlognvar(x(0.5*R*(R+1)+D+1)) ...
                + logpriormean(x(0.5*R*(R+1)+D+2));

logpost   = @(x) multitask_output_gp_loglik_vec( x, X, Y, all_difs2 ) + logprior(x); 

initial  = A(tril(ones(R))==1);
initial  = [initial;log(ls(:));log(nvar);mean];

samples  = slicesample(initial',nsamples,'logpdf',logpost,'burnin',20,'thin',5);

  A_samples  = {};
 ls_samples  = {};
nvar_samples = {};
Kinv_samples = {};
mean_samples = {};

difs2 = all_difs2(1:N,1:N,:);
difs2 = difs2 + tprod(difs2,[2,1,3],triu(ones(N)),[2,1]);

obj_mean_pred = zeros(N, R);
for ns=1:nsamples
    temp           = tril(ones(R));
    temp(temp==1)  = samples(ns,1:R*(R+1)/2);
    A_samples{ns}  = temp;
       
    ls_samples{ns} = exp(reshape(samples(ns,R*(R+1)/2+1:R*(R+1)/2+D),1,D));

    nvar_samples{ns} = exp(samples(R*(R+1)/2+D+1));    
    
    mean_samples{ns} = samples(R*(R+1)/2+D+2);    
    
    K       = nvar*eye(N);
    overls2 = 1./ls_samples{ns}./ls_samples{ns}; % 1 x D
    r2      = tprod(difs2,[1,2,-1],overls2',[-1]);  % N x N
    u       = sqrt(5*r2);   
    K       = K + (1+u+5/3*r2).*exp(-u);  %NxN
   
    Kinv_samples{ns} = solve_chol(chol(K),eye(N));
    
    K21 = K - nvar*eye(N);
    obj_mean_pred_ns = Kinv_samples{ns} * K21 * (Y - mean_samples{ns}) + mean_samples{ns};
    
    obj_mean_pred = obj_mean_pred + obj_mean_pred_ns / nsamples;
end


        





