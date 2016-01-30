% Y is NxR

function[ neglogEIV ] = neg_log_EIV_approx_all_cells_all_samples_noderiv( xstar, X, Y, A_samples, ls_samples, nvar_samples, Kinv_samples, mean_samples, cells )

f  = 0;
[N, D] = size(X);

nsamples = length(A_samples);
ncells = length(cells);

Yvec = Y';
Yvec = Yvec(:);

difs_star    = bsxfun(@times,ones(size(X,1),1),xstar) - X;     % N x D
difs2_star   = difs_star.*difs_star;

 logEIV_ns = zeros(nsamples, 1);
for ns=1:nsamples
    A       = A_samples{ns};
    ls      = ls_samples{ns};
    nvar    = nvar_samples{ns};
    K11inv  = Kinv_samples{ns};
    mean    = mean_samples{ns};

    logEIV_cells = zeros(ncells, 1);
    for s=1:ncells
        C = cells{s};
        Z   = 0.5*(  C(2,:)'-C(1,:)').^2;
        mu  = 0.5*(2*C(2,:)'+C(1,:)');
        tau = Z/9;

         g = log_EIV_approx_noderiv( xstar, X, Yvec, A, ls, nvar, Z, tau, mu, K11inv, mean, difs2_star );

         logEIV_cells(s)   =  g;
    end
    
     g = logsumexp(logEIV_cells);     % scalar
     p = exp(logEIV_cells - g);  % ncells x 1
    
     logEIV_ns(ns)    =  g;
end

 logEIV = logsumexp(logEIV_ns) - log(nsamples);    % scalar
      p = exp(logEIV_ns - logEIV + log(nsamples)); % nsamples x 1

  neglogEIV = -  logEIV;


