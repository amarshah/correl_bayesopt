
function [ X, Y, PF, PV ] = main( method, obj_file, nInitial, nSamples, nIter, savefile )

rng('shuffle')

addpath(genpath('./lightspeed'))
addpath(genpath('./TPROD'))
addpath('./utils')
addpath(genpath('./objectives'))
addpath(genpath('./methods'))

load(obj_file)

Xinit = repmat(xmin, nInitial, 1) + repmat((xmax - xmin), nInitial, 1) .* rand(nInitial, length(xmin));

if strcmp(method, 'correl')
    [ X,Y,PF,PV ] = runEIVcorrel( obj, noiseless_obj, lower_bound, upper_bound, Xinit, nSamples, nIter, xmin, xmax );    
elseif strcmp(method, 'indep')
    [ X,Y,PF,PV ] = runEIVindep( obj, noiseless_obj, lower_bound, upper_bound, Xinit, nSamples, nIter, xmin, xmax );        
elseif strcmp(method, 'multitask')
    [ X,Y,PF,PV ] = runEIVmultitask( obj, noiseless_obj, lower_bound, upper_bound, Xinit, nSamples, nIter, xmin, xmax );    
elseif strcmp(method, 'parego')
    [ X,Y,PF,PV ] = runParEGO( obj, noiseless_obj, lower_bound, upper_bound, Xinit, nSamples, nIter, xmin, xmax );    
elseif strcmp(method, 'random')
    [ X,Y,PF,PV ] = runRandom( obj, noiseless_obj, lower_bound, upper_bound, Xinit, nSamples, nIter, xmin, xmax );        
else
    disp('Unsupported method!')
    exit
end

save( savefile, 'X', 'Y', 'PF', 'PV' );

