% Y is NxR

function[ neglogEIV ] = neg_log_EIV_indep_all_cells_all_samples_noderiv( xstar, X, Y, amp_samples, ls_samples, nvar_samples, K22inv_samples, mean_samples, cells )

f  = 0;
[N, D] = size(X);

nsamples = length(amp_samples);

 logEIV_ns = zeros(nsamples, 1);

for ns=1:nsamples
    amp     = amp_samples{ns};
    ls      = ls_samples{ns};
    nvar    = nvar_samples{ns};
    K22invs = K22inv_samples{ns};
    mean    = mean_samples{ns};
    
    g = log_EIV_indep_all_cells_noderiv( xstar, X, Y, amp, ls, nvar, K22invs, mean, cells );

      logEIV_ns(ns)   =  g;
end

 logEIV = logsumexp(logEIV_ns) - log(nsamples);      % scalar

  neglogEIV = -  logEIV;

