
function [ X, Y, PF, PV ] = main( method, obj_file, nInitial, nSamples, nIter )

addpath('./utils')
addpath(genpath('./objectives'))

load(obj_file)

Xinit = repmat(xmin, nInitial, 1) + repmat((xmax - xmin), nInitial, 1) .* rand(nInitial, length(xmin));

if method=='correl'
    
    addpath(genpath('./methods/correl'))

    [ X,Y,PF,PV ] = runEIVcorrel( obj, noiseless_obj, lower_bound, upper_bound, Xinit, nSamples, nIter, xmin, xmax );
    
elseif method=='indep'
    
    addpath(genpath('./methods/indep'))
    
    [ X,Y,PF,PV ] = runEIVindep( obj, noiseless_obj, lower_bound, upper_bound, Xinit, nSamples, nIter, xmin, xmax );    
    
elseif method=='multitask'
    
    addpath(genpath('./methods/multitask'))

    [ X,Y,PF,PV ] = runEIVmultitask( obj, noiseless_obj, lower_bound, upper_bound, Xinit, nSamples, nIter, xmin, xmax );
    
elseif method=='parego'
    
    addpath(genpath('./methods/parego'))

    [ X,Y,PF,PV ] = runParEGO( obj, noiseless_obj, lower_bound, upper_bound, Xinit, nSamples, nIter, xmin, xmax );    
    
elseif method=='random'
    
    addpath('./methods/random')
    
    [ X,Y,PF,PV ] = runRandom( obj, noiseless_obj, lower_bound, upper_bound, Xinit, nSamples, nIter, xmin, xmax );    
    
else
    'Unsupported method!'
end



