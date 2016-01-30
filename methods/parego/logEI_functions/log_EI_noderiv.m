% computes the expected improvement at a new point

function[ f ] = log_EI_noderiv( xstar, X, Y, amp, ls, nvar, K22inv )

[ m_pred, var_pred ] = compute_pred_matrix_small_noderiv( xstar, X, Y, amp, ls, K22inv );

ybest = max(Y);

 std_pred = sqrt(var_pred);

 z = (m_pred - ybest) / std_pred;

 cdf_z = normcdf(z);
 pdf_z = normpdf(z);

 g = (m_pred - ybest)*cdf_z + std_pred*pdf_z;

 f = log(g);



