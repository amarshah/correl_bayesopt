
% from parego paper
% rescaled

function[ y ] = oka2( x )

x1 = x(:,1);
x2 = x(:,2);
x3 = x(:,3);

y1 = -1- x1 / 3.2;

y2 = (1.8 + (x1 + pi).^2/4/pi/pi ...
       - power(abs(x2-5*cos(x1)), 1/3) ...
       - power(abs(x3-5*sin(x1)), 1/3)) / 2.4;

y = [y1, y2]';    
    