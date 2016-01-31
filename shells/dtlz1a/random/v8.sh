#!/bin/sh
#
#$ -S /bin/bash
#
/usr/local/apps/matlab/matlabR2011b/bin/matlab - nodisplay -r "main( 'random', 'dtlz1a.mat', 5, 8, 100, '/bigscratch/as793/multi_objective_bayesopt/dtlz1a/random/v8.mat')"