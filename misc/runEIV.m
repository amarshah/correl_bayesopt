
% Y is NxR data matrix
% X is NxD
% kernel params: A is RxR, RD lengthscales, R noise params

% input: objectives, lower, upper, nSamples, nInitial, nIter, xmin, xmax, type

function[ X,Y,pareto_frontier,pareto_volume ] = runEIV( objectives, lower, upper, Xinit, nSamples, nIter, xmin, xmax, type ) 

D = length(xmin);
R = length(lower);

% sample random points 
nInitial = size(Xinit,1);
X = Xinit;
Y = zeros(nInitial,R);
for nI=1:nInitial
    y = objectives(X(nI,:));
    Y(nI,:) = y';
end

% compute pareto volume
pareto_volume    = [];
pareto_frontier  = {};
[ int_cells, pareto ]  = obs_to_integration_cells( Y, lower, upper );
pareto_volume(1)       = compute_pareto_volume( int_cells, lower, upper );
pareto_frontier{1}     = pareto;

% initialise samples
v = Y-repmat(mean(Y),nInitial,1);
v = sum(v.*v)/(nInitial-1);      %1xR
   A_samples{1} = diag(v);      %RxR
        temp    = 0.8*sqrt(xmax-xmin); %1xD
  ls_samples{1} = repmat(temp,R,1);  %RxD
nvar_samples{1} = 0.02*v;

for iter=1:nIter
    iter
    %   sample hyperparams                      
    [ A_samples, ls_samples, nvar_samples ] = sample_hypers( A_samples{end}, ls_samples{end}, nvar_samples{end}, X, Y, nSamples, type ); 

    %   optimize EIV                
    if strcmp(type,'correl')
        neg_acqui = @(x) neg_log_EIV_approx_all_cells_all_samples( x, X, Y, A_samples, ls_samples, nvar_samples, int_cells );
    elseif strcmp(type,'indep')
        neg_acqui = @(x) neg_log_EIV_indep_all_cells_all_samples( x, X, Y, A_samples, ls_samples, nvar_samples, int_cells );
    end        
    xnext = globalOptimization(neg_acqui, xmin, xmax);

    %   sample new point            
    ynext = objectives(xnext);
    X     = [X;xnext];
    Y     = [Y;ynext'];
    
    %   find pareto frontier, cells             
    [ int_cells, pareto ] = obs_to_integration_cells( Y, lower, upper );

    % compute volume of pareto frontier
    pareto_volume(iter+1)      = compute_pareto_volume( int_cells, lower, upper );
    pareto_frontier{iter+1}    = pareto;
end




