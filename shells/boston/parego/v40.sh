#!/bin/sh
#
#$ -S /bin/bash
#
/usr/local/apps/matlab/matlabR2011b/bin/matlab - nodisplay -r "main( 'parego', 'boston.mat', 5, 8, 100, '/bigscratch/as793/multi_objective_bayesopt/boston/parego/v40.mat')"