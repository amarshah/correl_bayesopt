% Y is NxR

function[neglogEIV] = neg_log_EIV_MCMC_all_cells_all_samples( xstar, X, Y, A_samples, ls_samples, nvar_samples, Kinv_samples, mean_samples, cells )

nsamples = length(A_samples);
[N, R] = size(Y);

Yvec = Y';
Yvec = Yvec(:);

difs_star    = bsxfun(@times,ones(size(X,1),1),xstar) - X;     % N x D
difs2_star   = difs_star.*difs_star;

EIV = 0;
for ns=1:nsamples
    A       = A_samples{ns};
    ls      = ls_samples{ns};
    nvar    = nvar_samples{ns};
    K11inv  = Kinv_samples{ns};
    mean    = mean_samples{ns};
    alldK12 = compute_derivs( xstar, X, A, ls );

    C = cells{1};
    Z   = 0.5*(  C(2,:)'-C(1,:)').^2;
    mu  = 0.5*(2*C(2,:)'+C(1,:)');
    tau = Z/9;
    
    [ ~, ~, ~, ~,...
        Kpredinv, logdetKpred, mpred,...
        ~, ~, ~, ~, ~, ~, ~ ] = compute_pred_matrices( xstar, X, A, ls, Yvec, tau, mu, K11inv, mean, difs2_star );
     
    for s=1:length(cells)
        C = cells{s};
        
        if R==3
            exp_vol = @(x,y,z) vol_pdf3(x,y,z,mpred,logdetKpred,Kpredinv,C(1,:)');
            g = integral3(exp_vol,C(1,1),C(2,1),C(1,2),C(2,2),C(1,3),C(2,3));
        elseif R==2
            exp_vol = @(x,y) vol_pdf2(x,y,mpred,logdetKpred,Kpredinv,C(1,:)');
            g = integral2(exp_vol,C(1,1),C(2,1),C(1,2),C(2,2));           
        end
        
        EIV = EIV + g / nsamples;        
    end
end

neglogEIV = -log(EIV);

