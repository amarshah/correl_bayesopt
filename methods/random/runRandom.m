%%%%%
% Models scalar combination of objectives as a GP
%%%%%

% Y is NxR data matrix
% X is NxD

% input: objectives, lower, upper, nSamples, nInitial, nIter, xmin, xmax

function[ X,Y,pareto_frontier,pareto_volume ] = runRandom( objectives, noiseless_objectives, lower, upper, Xinit, nSamples, nIter, xmin, xmax ) 

D = length(xmin);
R = length(lower);

% sample random points 
N = size(Xinit,1);
X = Xinit;
Y = zeros(N,R);
F = zeros(N,R);
for nI=1:N
    y = objectives(X(nI,:));
    Y(nI,:) = y';
%    f = noiseless_objectives(X(nI,:));
    F = Y;
%    F(nI,:) = f';
end

pareto_volume    = [];
pareto_frontier  = {};

for iter=1:nIter
    iter
    
    %   find pareto volume based on noiseless objectives at probed locations
    [ F_cells, pareto ] = obs_to_integration_cells( F, lower, upper );    
    pareto_volume(iter)   = compute_pareto_volume( F_cells, lower, upper );
    pareto_frontier{iter} = pareto;

    % choose new point 
    xnext = rand(1,D).*(xmax - xmin) + xmin;
    
    %   sample new point            
    ynext = objectives(xnext);
%    fnext = noiseless_objectives(xnext);
    X     = [X;xnext];
    Y     = [Y;ynext'];
    F     = Y;
%    F     = [F;fnext'];
    N     = N+1;
    
end

%   find pareto volume based on noiseless objectives at probed locations
[ F_cells, pareto ] = obs_to_integration_cells( F, lower, upper );    
pareto_volume(nIter+1)   = compute_pareto_volume( F_cells, lower, upper );
pareto_frontier{nIter+1} = pareto;


