#!/bin/sh
#
#$ -S /bin/bash
#
/usr/local/apps/matlab/matlabR2013a/bin/matlab - nodisplay -r "main( 'correl', 'rocket.mat', 5, 8, 100, '/bigscratch/as793/multi_objective_bayesopt/rocket/correl/v11.mat')"