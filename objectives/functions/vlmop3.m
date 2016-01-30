% from ParEgo paper
% rescaled to have similar scaled outputs

function[ y ] = vlmop3(x)

x1 = x(:,1);
x2 = x(:,2);

x2plusy2 = x1.^2 + x2.^2;

y1 = -3-(0.5*(x2plusy2) + sin(x2plusy2))/4.1;

y2 = 0.5-((3*x1-2*x2+4).^2/8 + ...
       (x1-x2+1).^2/27 -24)/50; 

y3 = 3.3*(0.05 - 1./(x2plusy2+1) + 1.1*exp(-x2plusy2));

y = [y1,y2,y3]';


