#!/bin/sh
#
#$ -S /bin/bash
#
/usr/local/apps/matlab/matlabR2013a/bin/matlab - nodisplay -r "main( 'random', 'rocket.mat', 5, 8, 100, '/bigscratch/as793/multi_objective_bayesopt/rocket/random/v38.mat')"