
function[ f ] = neg_log_EI_all_samples_noderiv( xstar, X, Y, amp_samples, ls_samples, nvar_samples, K22inv_samples)

f  = 0;

nsamples = length(amp_samples);

for ns=1:nsamples
    amp     = amp_samples{ns};
    ls      = ls_samples{ns};
    nvar    = nvar_samples{ns};
    K22inv  = K22inv_samples{ns};
    
    g = log_EI_noderiv( xstar, X, Y, amp, ls, nvar, K22inv );

    f =  f -  g;    
end

