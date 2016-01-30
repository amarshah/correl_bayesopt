% rocket_der.m: returns the derivatives for the rocket equations
function derivs = rocket_der( t, w, flag,alpha,ux,uy,M,tfmax)
% w(1):x, w(2):vx, w(3):y, w(4):vy, w(5):m
% fuel runs out after tfmax, so use a step function to stop burn rate alpha.
wr=sqrt(w(1).^2+w(3).^2);
tmp=alpha*stepf(t,tfmax,1);
derivs = [w(2);ux*tmp./w(5)-M*w(1)./wr.^3;w(4);...
               uy*tmp./w(5)-M*w(3)./wr.^3;-tmp;];
