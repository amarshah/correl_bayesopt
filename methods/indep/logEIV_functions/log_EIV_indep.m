% xstar   D  x 1
% X       N  x D
% Y       N  x R 
% amp     1  x R
% ls      RxD lengthscales
% nvar    Rx1 noise vars
% mean    
% lower   1  x R
% upper   1  x R

%  logEIV     scalar
% dlogEIV     Dx1
function[ logEIV, dlogEIV ] = log_EIV_indep( xstar, X, Y, amp, ls, nvar, K22inv, mean, lower, upper )

[N,D] = size(X);
 R    = length(amp);

 logEIV = 0;
dlogEIV = zeros(D,1);
for r=1:R
    Yr = Y(:,r);
    [ m_pred, v_pred, dm_pred, dv_pred ] = compute_pred_matrix_small( xstar, X, Yr, amp(r), ls(r,:), K22inv{r} );
    
    lr = lower(r);
    ur = upper(r);

     a = (ur-m_pred)/sqrt(v_pred); 
     b = (lr-m_pred)/sqrt(v_pred); 
    da = -dm_pred/sqrt(v_pred) - 0.5*a/v_pred*dv_pred; 
    db = -dm_pred/sqrt(v_pred) - 0.5*b/v_pred*dv_pred; 

    cdf_a = normcdf(a);
    cdf_b = normcdf(b);
    pdf_a = normpdf(a);
    pdf_b = normpdf(b);
    
     temp     = sqrt(v_pred)*(pdf_b-pdf_a) + (m_pred-lr)*(cdf_a-cdf_b);
    dtemp     = -0.5*temp/v_pred*dv_pred + sqrt(v_pred)*( da*pdf_a-db*pdf_b ) + ...
                 dm_pred*(cdf_a-cdf_b) + (m_pred-lr)*(da*pdf_a-db*pdf_b);
     logEIV   = logEIV + log(temp);
    dlogEIV   = dlogEIV + dtemp/temp;
end




