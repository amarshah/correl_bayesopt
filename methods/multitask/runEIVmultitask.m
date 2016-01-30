%%%%%%%%%%%%%%%
% Assumes a kronecker product kernel accross tasks
%%%%%%%%%%%%%%%

% Y is NxR data matrix
% X is NxD
% kernel params: A is RxR, RD lengthscales, R noise params

% input: objectives, lower, upper, nSamples, nInitial, nIter, xmin, xmax

function[ X,Y,pareto_frontier,pareto_volume ] = runEIVmultitask( objectives, noiseless_objectives, lower, upper, Xinit, nSamples, nIter, xmin, xmax ) 

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

% initialise samples
m = sum(sum(Y))/N/R;
v = Y-m*ones(size(Y));
v = sum(sum(v.*v))/(N*R-1);      %scalar
  ls_samples{1} = 0.5*sqrt(xmax-xmin); %1xD
nvar_samples{1} = 0.04*v;
A_samples{1} = .85*v*eye(R);
mean_samples{1} = m;

pareto_volume    = [];
pareto_frontier  = {};

for iter=1:nIter
    iter
    %   sample hyperparams                      
    [ A_samples, ls_samples, nvar_samples, Kinv_samples, mean_samples, obj_means ] = sample_hypers_multitask( A_samples{end}, ls_samples{end}, nvar_samples{end}, mean_samples{end}, X, difX2, Y, nSamples, lower, upper );     
    
    %   find integration cells based on prediction of noiseless objectives             
    [ int_cells, ~ ] = obs_to_integration_cells( Y, lower, upper );

    %   find pareto volume based on noiseless objectives at probed locations
    [ F_cells, pareto ] = obs_to_integration_cells( F, lower, upper );    
    pareto_volume(iter)   = compute_pareto_volume( F_cells, lower, upper );
    pareto_frontier{iter} = pareto;
    
    %   optimize EIV                
    neg_acqui_noderiv = @(x) neg_log_EIV_multitask_all_cells_all_samples_noderiv( x, X, Y, A_samples, ls_samples, nvar_samples, Kinv_samples, mean_samples, int_cells );
    neg_acqui = @(x) neg_log_EIV_multitask_all_cells_all_samples( x, X, Y, A_samples, ls_samples, nvar_samples, Kinv_samples, mean_samples, int_cells );
    xnext = globalOptimization(neg_acqui_noderiv, neg_acqui, xmin, xmax);

    %   sample new point            
    ynext = objectives(xnext);
    X     = [X;xnext];
    Y     = [Y;ynext'];
    F     = Y;
    N     = N+1;
    
    % update difX2
    temp = bsxfun(@times,ones(N-1,1),xnext) - X(1:N-1,:);     % N-1 x D    
    difX2(1:N-1,N,:) = reshape(temp.*temp,N-1,1,D);        
end

%   find pareto volume based on noiseless objectives at probed locations
[ F_cells, pareto ] = obs_to_integration_cells( F, lower, upper );    
pareto_volume(nIter+1)   = compute_pareto_volume( F_cells, lower, upper );
pareto_frontier{nIter+1} = pareto;

