
% mpred        Rx1
% Kinv         RxR
% cell_floor   Rx1

function[ f ] = vol_pdf3( x, y, z, mpred, logdetK, Kinv, cell_floor )

R = length(cell_floor);
[A, B] = size(x);

lognormconst = - 0.5*R*log(2*pi) - 0.5*logdetK;
f = zeros(A,B);
for a=1:A
    for b=1:B
        X = [x(a,b), y(a,b), z(a,b)]';
        
        dif = X - mpred;
        logpdf = lognormconst - 0.5*dif'*Kinv*dif;  % Mx1
        pdf = exp(logpdf);
        vol = prod(X-cell_floor);
        
        f(a,b) = pdf*vol;
    end
end
