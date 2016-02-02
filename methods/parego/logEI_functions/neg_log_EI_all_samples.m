
function[ f, df ] = neg_log_EI_all_samples( xstar, X, Y, amp_samples, ls_samples, nvar_samples, K22inv_samples)

f  = 0;
df = zeros(length(xstar),1);

nsamples = length(amp_samples);

for ns=1:nsamples
    amp     = amp_samples{ns};
    ls      = ls_samples{ns};
    nvar    = nvar_samples{ns};
    K22inv  = K22inv_samples{ns};
    
    [ g, dg ] = log_EI( xstar, X, Y, amp, ls, nvar, K22inv );
    
    if isreal(g)==0
       g; 
    end

    f =  f -  g;
    df = df - dg;
end

