
function[ mpred ] = gp_minushalf(xstar)

load minushalf

N = size(X,1);
temp = bsxfun(@times,ones(N,1),xstar) - X; 
difs2_star = temp.*temp;

tau = 1; mu = 0; 
[ ~,~,~,~,~,~,mpred,~,~,~,~,~,~,~] = compute_pred_matrices( xstar, X, A, lss, y, tau, mu, Kinv, meam, difs2_star );


