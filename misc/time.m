
r = randn(700, 700);

u = r'*r;

A = chol(u);
B = cholproj(u);
