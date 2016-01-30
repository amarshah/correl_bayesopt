
% takes inputs X (NxD) and multiobj outputs (NxR) and creates a random
% linear combination of objectives based on eq 2 of ParEgo paper

function[ newY ] = parego_scalarize(X, Y, lower, upper)

[N, R] = size(Y);

% first map observations to unit scale output 
scaledY = (Y - bsxfun(@times, ones(N,1), lower)) ./ bsxfun(@times, ones(N,1), upper-lower);

% sample weights
lambda = gamrnd(1,1,1,R);
lambda = lambda ./ sum(lambda);

% compute augmented Tchebycheff function
lambda_scaledY = bsxfun(@times, ones(N,1), lambda) .* scaledY;   % NxR

newY = min(scaledY')' + 0.05 * sum(scaledY')';   % Nx1



