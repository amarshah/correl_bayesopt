
l = [-10,-10];
u = [40,30];

m = [7,9];
K = [2,1 ;...
     1,3 ];

fun2 = @(x,y) integrand(x,y,l,m,K);

ev  = integral2(fun2,l(1),u(1),l(2),u(2))

Zs    = (u-l).^2/2;
means = (2*u+l)/3;
vars  = (u-l).^2/18;

Kinv  = inv(K);
Sigma = inv(Kinv + diag(1./vars));
Mu    = (m*Kinv + means./vars)*Sigma;

logZ  = -0.5*(m*Kinv*m'+log(det(K)))+...
          sum(log(Zs)-0.5*(means.^2./vars+log(vars)+log(2*pi)))+...
          0.5*(Mu*inv(Sigma)*Mu'+log(det(Sigma)));
approx = exp(logZ)
      
 
 