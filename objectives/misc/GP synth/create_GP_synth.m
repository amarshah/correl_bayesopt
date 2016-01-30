
   A = [1, 0; 0, 1 ]; 
meam = [-2, 2];
  
 lss = [0.7, 0.4, 1;
        0.34, 0.9, 0.5];

nvar = [0.01,0.01];

   R = 2;
   
xmin = [0,0,0];
xmax = [1,1,1];
D = length(xmin);


N = 400;
X = rand(N,D);
difs2 = zeros(N,N,D);
for n=1:N-1
    temp = bsxfun(@times,ones(N-n,1),X(n,:)) - X(n+1:N,:); 
    difs2(n,n+1:N,:) = reshape(temp.*temp,1,N-n,D);
end


K = compute_kernel( X, difs2, A, lss, nvar );
L = chol(K);
y = L'*randn(R*N, 1) + repmat(meam', N, 1);

C    = chol(K);
Kinv = solve_chol(C,eye(N*R));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xstar = rand(1,D);

temp = bsxfun(@times,ones(N,1),xstar) - X; 
difs2_star = temp.*temp;

tau = 1; mu = 0; 
[ ~,~,~,~,~,~,mpred,~,~,~,~,~,~,~] = compute_pred_matrices( xstar, X, A, lss, y, tau, mu, Kinv, meam, difs2_star );

%save('zero.mat','X','A','lss','y','Kinv','nvar', 'meam');


