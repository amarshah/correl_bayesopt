#!/bin/sh
#
#$ -S /bin/bash
#
/usr/local/apps/matlab/matlabR2011b/bin/matlab - nodisplay -r "main( 'random', 'oka2.mat', 5, 8, 100, '/bigscratch/as793/multi_objective_bayesopt/oka2/random/v47.mat')"