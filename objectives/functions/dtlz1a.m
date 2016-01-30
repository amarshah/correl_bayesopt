
% from parego paper
% rescaled

function[ y ] = dtz1a( x )

x1 = x(:,1);

g = 500;

for i=2:6
    g = g + (x(:,i)-0.5).^2 - cos(2*pi*(x(:,i)-0.5)); 
end

y1 = -0.5.*x1.*(1+g)/128 - 1.01;

y2 = -0.5.*(1-x1).*(1+g)/128 + 2.99;

y = [y1, y2]';
