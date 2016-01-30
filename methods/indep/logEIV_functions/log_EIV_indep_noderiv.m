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
function[ logEIV ] = log_EIV_indep_noderiv( xstar, X, Y, amp, ls, nvar, K22inv, mean, lower, upper )

[N,D] = size(X);
 R    = length(amp);

 logEIV = 0;
for r=1:R
    Yr = Y(:,r);
    [ m_pred, v_pred ] = compute_pred_matrix_small_noderiv( xstar, X, Yr, amp(r), ls(r,:), K22inv{r} );
    
    lr = lower(r);
    ur = upper(r);

     a = (ur-m_pred)/sqrt(v_pred); 
     b = (lr-m_pred)/sqrt(v_pred); 

    cdf_a = normcdf(a);
    cdf_b = normcdf(b);
    pdf_a = normpdf(a);
    pdf_b = normpdf(b);
    
     temp     = sqrt(v_pred)*(pdf_b-pdf_a) + (m_pred-lr)*(cdf_a-cdf_b);
     logEIV   = logEIV + log(temp);
end




