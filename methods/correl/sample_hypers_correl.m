% A       R  x R
% ls      RxD lengthscales
% nvar    Rx1 noise vars
% Y       NxR  

function[ A_samples, ls_samples, nvar_samples, Kinv_samples, mean_samples, obj_mean_pred ] = sample_hypers_correl( A, ls, nvar, mean, X, all_difs2, Y, nsamples, lower, upper )
  
[N,D] = size(X);
 R    = size(A,1);

% hyper priors: 
% A         elements are N( mA , vA  )
% logls     elements are N( mls , vls )
% lognvar   elements are N( mvar , vvar   ) 
% mean      elements are N( 0, meanvar )

vA   = 0.5;    mA   = 0;  
vls  = 0.3;   mls  = -0.6;
mvar = log((upper-lower).^2/40/50); vvar = abs(mvar)/2; 
meanmean = (upper+lower)/2;  meanvar = (upper-lower).^2/40;

logpriorA       = @(x) -0.5*sum((x-mA).*(x-mA))/vA;
logpriorlogls   = @(x) -0.5*sum((x-mls).*(x-mls))/vls - sum(x);   % '-x' is for the log change of var
logpriorlognvar = @(x) -0.5*sum((x(:)-mvar(:)).*(x(:)-mvar(:))./vvar(:)) - sum(x);
logpriormean    = @(x) -0.5*sum((x(:)-meanmean(:)).*(x(:)-meanmean(:))./meanvar(:));

%sum(log(1e-300+double(lower'<=x & x<=upper')));

% first R*(R+1)/2 elt are for A,
% next  R*D elts are for logls,
% next  R elts are for lognvar,
% next  R elts are for mean,
logprior  = @(x) logpriorA(x(1:0.5*R*(R+1))) ...
                + logpriorlogls(x(0.5*R*(R+1)+1:0.5*R*(R+1)+R*D)) ...
                + logpriorlognvar(x(0.5*R*(R+1)+R*D+1:0.5*R*(R+1)+R*D+R)) ...
                + logpriormean(x(0.5*R*(R+1)+R*D+R+1:end));
            
logpost   = @(x) correl_output_gp_loglik_vec( x, X, Y, all_difs2 ) + logprior(x); 

initial  = [A(tril(ones(R))==1);log(ls(:));log(nvar(:));mean(:)];

samples  = slicesample(initial',nsamples,'logpdf',logpost,'burnin',20,'thin',5);

  A_samples  = {};
 ls_samples  = {};
nvar_samples = {};
Kinv_samples = {};
mean_samples = {};

Yvec   = Y';
Yvec   = Yvec(:);
obj_mean_pred = zeros(N, R);

for ns=1:nsamples
    temp           = tril(ones(R));
    temp(temp==1)  = samples(ns,1:R*(R+1)/2);
    A_samples{ns}  = temp;
       
    ls_samples{ns} = exp(reshape(samples(ns,R*(R+1)/2+1:R*(R+1)/2+R*D),R,D));

    nvar_samples{ns} = max(exp(samples(ns,R*(R+1)/2+R*D+1:R*(R+1)/2+R*D+R)),1e-6);    
    
    mean_samples{ns} = samples(ns,R*(R+1)/2+R*D+R+1:end);
    
    K = compute_kernel( X, all_difs2, A_samples{ns}, ls_samples{ns}, nvar_samples{ns} );
    Kinv_samples{ns} = solve_chol(chol(K),eye(N*R));
    
    K21 = K - diag(repmat(nvar_samples{ns}, 1, N));
    
    mean_rep = repmat(mean_samples{ns}', N, 1);
    K22invK21 = Kinv_samples{ns} * K21;    % N*R x N*R
    obj_mean_pred_ns = K22invK21 * (Yvec - mean_rep) + mean_rep;  % N*R x 1 
    
    obj_mean_pred_ns = reshape(obj_mean_pred_ns, R, N)';   % N x R
    
    obj_mean_pred = obj_mean_pred + obj_mean_pred_ns / nsamples; 

end

