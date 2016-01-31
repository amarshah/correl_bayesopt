% Y is NxR

function[f,df] = log_EIV_approx_all_cells( xstar, X, Y, A, ls, nvar, cells )

f  = 0;
df = zeros(length(xstar),1);

Yvec = Y';
Yvec = Yvec(:);

for s=1:length(cells)
    C = cells{s};
    Z   = 0.5*(  C(2,:)-C(1,:)).^2;
    mu  = 0.5*(2*C(2,:)+C(1,:));
    tau = Z/9;

    [ g, dg ] = log_EIV_approx( xstar, X, Yvec, A, ls, nvar, Z, tau, mu );
    
     f =  f +  g;
    df = df + dg;
end



