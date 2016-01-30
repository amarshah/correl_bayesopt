% amp     R  x 1
% ls      RxD lengthscales
% nvar    Rx1 noise vars
% Y       NxR  

function[ amp_samples, ls_samples, nvar_samples, K22inv_samples, mean_samples, obj_mean_pred ] = sample_hypers_indep( amp, ls, nvar, mean, X, all_difs2, Y, nsamples, lower, upper )
  
[N,D] = size(X);
 R    = length(amp);
 
% hyper priors: 
% logamp    elements are N( mamp , vamp  )
% logls     elements are N( mls , vls )
% lognvar   elements are N( mvar , vvar   ) 
% mean      elements are N( meanmean, meanvar )

vamp = 0.75;   mamp = 0.2;  
vls  = 0.5;   mls  = -0.6;
%vvar = 0.5;   mvar = -2.75; 
mvar = log((upper-lower).^2/40/50); vvar = abs(mvar)/2; 
meanmean = (upper+lower)/2;  meanvar = (upper-lower).^2/40;

logpriorlogamp  = @(x) -0.5*sum((x-mamp).*(x-mamp))/vamp - sum(x); % '-x' is for the log change of var
logpriorlogls   = @(x) -0.5*sum((x-mls).*(x-mls))/vls - sum(x);
logpriorlognvar = @(x) -0.5*sum((x(:)-mvar(:)).*(x(:)-mvar(:))./vvar(:)) - sum(x);
%logpriorlognvar = @(x) -0.5*sum((x-mvar).*(x-mvar))/vvar - sum(x);
logpriormean    = @(x) -0.5*sum((x-meanmean).*(x-meanmean)./meanvar);

% first R elts are for logamp,
% next  R*D elts are for logls,
% next  R elts are for lognvar,
% next  R elts are for mean,
logprior  = @(x) logpriorlogamp(x(1:R)) ...
                + logpriorlogls(x(R+1:R+R*D)) ...
                + logpriorlognvar(x(R+R*D+1:R+R*D+R)) ...
                + logpriormean(x(2*R+R*D+1:end));

logpost   = @(x) indep_output_gp_loglik_vec( x, X, Y, all_difs2 ) + logprior(x); 

initial  = [log(amp(:));log(ls(:));log(nvar(:));mean(:)];

samples  = slicesample(initial',nsamples,'logpdf',logpost,'burnin',20,'thin',5);

amp_samples  = {};
 ls_samples  = {};
nvar_samples = {};
mean_samples = {};
K22inv_samples = {};

obj_mean_pred = zeros(N, R);

for ns=1:nsamples
    amp_samples{ns} = exp(samples(ns,1:R));
       
    ls_samples{ns}  = exp(reshape(samples(ns,R+1:R+R*D),R,D));

    nvar_samples{ns} = exp(samples(ns,R+R*D+1:2*R+R*D));    
    
    mean_samples{ns} = samples(ns,2*R+R*D+1:end);    
    
    K22invs = zeros(N,N,R);
    K21s    = zeros(N,N,R);
    for r=1:R
        amp  =  amp_samples{ns}(r);
        ls   =   ls_samples{ns}(r,:);
        nvar = nvar_samples{ns}(r);

        K = (amp+nvar)*eye(N);
        
        overls2 = 1./ls./ls;   % 1 x D
        for n=1:N-1
            difs2         = reshape(all_difs2(n,n+1:N,:),N-n,D);
            r2            = difs2*overls2';  %N-n x1
            u             = sqrt(5*r2);
            K(n+1:N,n)    = amp*(1+u+5/3*r2).*exp(-u);  
            K(n,n+1:N)    = K(n+1:N,n)'; 
        end
        
        K21s(:,:,r) = K - nvar*eye(N);
        
        C          = chol(K);
        K22invs(:,:,r) = solve_chol(C,eye(N));
    end
    
    K22inv_samples{ns} = K22invs;
    
    % get contribution to the mean prediction
    mean_tile = bsxfun(@times, ones(N,1), mean_samples{ns});   % N x R
    K22invK21 = tprod(K22invs, [1, -1, 3], K21s, [-1, 2, 3]);    % N x N x R
    obj_mean_pred_ns = tprod(K22invK21, [1, -1, 2], Y - mean_tile, [-1, 2]) + mean_tile;  % N x R 
    
    obj_mean_pred = obj_mean_pred + obj_mean_pred_ns / nsamples; 
end

        
end




