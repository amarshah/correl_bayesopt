
function[ y ] = multi_test( x )

y = zeros(3,1);
y(1) = x(1)*x(2) - cos(x(4)*x(3))+1;
y(2) = x(5)*sqrt(x(3)) + 2*sin(x(4))+2;
y(3) = x(4)*sqrt(x(1)) + 2*sin(x(3))+2;

