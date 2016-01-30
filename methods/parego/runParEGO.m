%%%%%
% Models scalar combination of objectives as a GP
%%%%%

% Y is NxR data matrix
% X is NxD

% input: objectives, lower, upper, nSamples, nInitial, nIter, xmin, xmax

function[ X,Y,pareto_frontier,pareto_volume ] = runParEGO( objectives, noiseless_objectives, lower, upper, Xinit, nSamples, nIter, xmin, xmax ) 

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

% compute difX2
difX2 = zeros(N+nIter,N+nIter,D);
for n=1:N
    temp = bsxfun(@times,ones(N-n,1),X(n,:)) - X(n+1:N,:);     % N-n x D    
    difX2(n,n+1:N,:) = reshape(temp.*temp,1,N-n,D);
end 

pareto_volume    = [];
pareto_frontier  = {};

for iter=1:nIter
    iter
    %   sample outputs of augmented Tchebycheff function (eq 2 of ParEgo paper)
    newY = parego_scalarize(X, Y, lower, upper);

    %   sample hyperparams                   
    [ amp_samples, ls_samples, nvar_samples, Kinv_samples, mean_samples ] = sample_hypers_parego( X, difX2, newY, nSamples, lower, upper );

    %   find pareto volume based on noiseless objectives at probed locations
    [ F_cells, pareto ] = obs_to_integration_cells( F, lower, upper );    
    pareto_volume(iter)   = compute_pareto_volume( F_cells, lower, upper );
    pareto_frontier{iter} = pareto;
    
    %   optimize EI of scalarized objective
    neg_acqui_noderiv = @(x) neg_log_EI_all_samples_noderiv( x, X, newY, amp_samples, ls_samples, nvar_samples, Kinv_samples);
    neg_acqui = @(x) neg_log_EI_all_samples( x, X, newY, amp_samples, ls_samples, nvar_samples, Kinv_samples);
    xnext = globalOptimization(neg_acqui_noderiv, neg_acqui, xmin, xmax);
    
    %   sample new point            
    ynext = objectives(xnext);
%    fnext = noiseless_objectives(xnext);
    X     = [X;xnext];
    Y     = [Y;ynext'];
    F     = Y;
%    F     = [F;fnext'];
    N     = N+1;
    
    %   update difX2
    temp = bsxfun(@times,ones(N-1,1),xnext) - X(1:N-1,:);     % N-1 x D    
    difX2(1:N-1,N,:) = reshape(temp.*temp,N-1,1,D);            
end

%   find pareto volume based on noiseless objectives at probed locations
[ F_cells, pareto ] = obs_to_integration_cells( F, lower, upper );    
pareto_volume(nIter+1)   = compute_pareto_volume( F_cells, lower, upper );
pareto_frontier{nIter+1} = pareto;


