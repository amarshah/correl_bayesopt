% computes the expected improvement at a new point

function[ f, df ] = log_EI( xstar, X, Y, amp, ls, nvar, K22inv )

[ m_pred, var_pred, dm_pred, dvar_pred ] = compute_pred_matrix_small( xstar, X, Y, amp, ls, K22inv );

ybest = max(Y);

 std_pred = sqrt(var_pred);
dstd_pred = 0.5*std_pred*dvar_pred / var_pred;

 z = (m_pred - ybest) / std_pred;
dz = dm_pred/std_pred - z*dstd_pred/std_pred;

 cdf_z = normcdf(z);
 pdf_z = normpdf(z);

dcdf_z = pdf_z * dz;
dpdf_z = -z * dcdf_z;

 g = (m_pred - ybest)*cdf_z + std_pred*pdf_z;
dg = (m_pred - ybest)*dcdf_z + dm_pred*cdf_z ...
        + dstd_pred*pdf_z + std_pred*dpdf_z;

if g<0
     f = 0; 
    df = zeros(size(dg)); 
else    
     f = log(g);
    df = dg / (g+1e-10);
end

