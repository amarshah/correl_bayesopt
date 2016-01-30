
% from parego paper
% rescaled

function[ y ] = dtz2a( x )

x1 = x(:,1);
x2 = x(:,2);
d = size(x,2);

g = 0;

for i=3:d
    g = g + (x(:,i)-0.5).^2; 
end

y1 = 1.04-(1+g).*cos(x1*pi/2).*cos(x2*pi/2);

y2 = 1.04-(1+g).*cos(x1*pi/2).*cos(x2*pi/2);

y3 = 1.1-(1+g).*sin(x1*pi/2);

y = [y1, y2, y3]'/1.16;

y(1) = y(1) + 2;
y(2) = y(2) - 2;
y(3) = y(3) / 4;

