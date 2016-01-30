
function[ y ] = multi_test_noisy( x )

y = multi_test( x );

y(1) = y(1) + 0.10 * randn;
y(2) = y(2) + 0.06 * randn;
y(3) = y(3) + 0.14 * randn;
