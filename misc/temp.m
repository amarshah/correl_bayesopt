
function[ y dy ] = temp(x)

    difs = x - [0,0,0]';     % Dx1

    ls      = [1.1,0.8,1.6]';
    temp    = 1./ls./ls;
    r2      = (difs.*difs)'*temp;  % scalar
    dr2     = 2*difs./ls./ls;           % Dx1

    dy = -5/6*dr2.*(1+sqrt(5*r2)).*exp(-sqrt(5*r2)); 
     y = (1+sqrt(5*r2)+5/3*r2).*exp(-sqrt(5*r2)); 

