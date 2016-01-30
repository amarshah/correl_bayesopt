clear

objectives = @(x) multi_test_noisy(x);
noiseless_objectives = @(x) multi_test(x);

xmin = [0,0,0,0,0];
xmax = [1,1,1,1,1];

nSamples = 8;
nInitial = 5;
nIter    = 3;

Xinit = repmat(xmin, nInitial, 1) + repmat((xmax - xmin), nInitial, 1) .* rand(nInitial, length(xmin));

lower = [-1,-1,-1];
upper = [4,6,6];

[ X2,Y2,pareto_frontier2,pareto_volume2 ] = runEIVcorrel( objectives, noiseless_objectives, lower, upper, Xinit, nSamples, nIter, xmin, xmax); 
[ X1,Y1,pareto_frontier1,pareto_volume1 ] = runEIVindep(  objectives, noiseless_objectives, lower, upper, Xinit, nSamples, nIter, xmin, xmax); 
[ X4,Y4,pareto_frontier4,pareto_volume4 ] = runParEGO(  objectives, noiseless_objectives, lower, upper, Xinit, nSamples, nIter, xmin, xmax); 
[ X3,Y3,pareto_frontier3,pareto_volume3 ] = runEIVmultitask( objectives, noiseless_objectives, lower, upper, Xinit, nSamples, nIter, xmin, xmax); 

%%%%%%%%%%
% still getting negative pareto volumes for some reason 

