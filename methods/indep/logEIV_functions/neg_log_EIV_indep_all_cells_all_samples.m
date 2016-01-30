% Y is NxR

function[ neglogEIV , dneglogEIV ] = neg_log_EIV_indep_all_cells_all_samples( xstar, X, Y, amp_samples, ls_samples, nvar_samples, K22inv_samples, mean_samples, cells )

f  = 0;
df = zeros(length(xstar),1);
[N, D] = size(X);

nsamples = length(amp_samples);

 logEIV_ns = zeros(nsamples, 1);
dlogEIV_ns = zeros(nsamples, D);

for ns=1:nsamples
    amp     = amp_samples{ns};
    ls      = ls_samples{ns};
    nvar    = nvar_samples{ns};
    K22invs = K22inv_samples{ns};
    mean    = mean_samples{ns};
    
    [ g, dg ] = log_EIV_indep_all_cells( xstar, X, Y, amp, ls, nvar, K22invs, mean, cells );

      logEIV_ns(ns)   =  g;
     dlogEIV_ns(ns,:) = dg;
end

 logEIV = logsumexp(logEIV_ns) - log(nsamples);      % scalar
      p = exp(logEIV_ns - logEIV + log(nsamples));   % nsamples x 1
dlogEIV = sum(bsxfun(@times, p, dlogEIV_ns));        % 1 x D  

  neglogEIV = -  logEIV;
 dneglogEIV = - dlogEIV;

