% sort derivs and summing over cells

% Y is NxR

function[ neglogEIV , dneglogEIV ] = neg_log_EIV_multitask_all_cells_all_samples( xstar, X, Y, A_samples, ls_samples, nvar_samples, Kinv_samples, mean_samples, cells )

f  = 0;
df = zeros(length(xstar),1);
[N, D] = size(X);

nsamples = length(A_samples);
ncells = length(cells);

% Yvec = Y';
% Yvec = Yvec(:);

difs_star    = bsxfun(@times,ones(size(X,1),1),xstar) - X;     % N x D
difs2_star   = difs_star.*difs_star;

 logEIV_ns = zeros(nsamples, 1);
dlogEIV_ns = zeros(nsamples, D);

for ns=1:nsamples
    A       = A_samples{ns};
    ls      = ls_samples{ns};
    nvar    = nvar_samples{ns};
    K11inv  = Kinv_samples{ns};
    mean    = mean_samples{ns};
    G       = A*A';
    Ginv    = solve_chol(chol(G),eye(size(G,1)));
    
    overls2 = 1./ls./ls;
         r2 = difs2_star*overls2';   %Nx1
          u = sqrt(5*r2);     
        K12 = (1+u+5/3*r2).*exp(-u); % Nx1

        dr2 = 2*bsxfun(@times,difs2_star,overls2);    %NxD
          v = (1+u).*exp(-u);
       dK12 = -5/6*bsxfun(@times,v,dr2);  %NxD
    
     logEIV_cells = zeros(ncells, 1);
    dlogEIV_cells = zeros(ncells, D);
    for s=1:ncells
        C = cells{s};
        Z   = 0.5*(  C(2,:)'-C(1,:)').^2;
        mu  = 0.5*(2*C(2,:)'+C(1,:)');
        tau = Z/9;

        [ g, dg ] = log_EIV_multitask( xstar, Y, Z, tau, mu, G, Ginv, K11inv, mean, K12, dK12 );

         logEIV_cells(s)   =  g;
        dlogEIV_cells(s,:) = dg; 
    end

     g = logsumexp(logEIV_cells);                % scalar
     p = exp(logEIV_cells - g);                  % ncells x 1
    dg = sum(bsxfun(@times, p, dlogEIV_cells));  % 1 x D  
    
     logEIV_ns(ns)    =  g;
    dlogEIV_ns(ns, :) = dg;

end

 logEIV = logsumexp(logEIV_ns) - log(nsamples);      % scalar
      p = exp(logEIV_ns - logEIV + log(nsamples));   % nsamples x 1
dlogEIV = sum(bsxfun(@times, p, dlogEIV_ns));        % 1 x D  

  neglogEIV = -  logEIV;
 dneglogEIV = - dlogEIV;


