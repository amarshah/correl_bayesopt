#!/bin/sh
#
#$ -S /bin/bash
#
/usr/local/apps/matlab/matlabR2011b/bin/matlab - nodisplay -r "main( 'indep', 'gp_zero.mat', 5, 8, 100, '/bigscratch/as793/multi_objective_bayesopt/gp_zero/indep/v28.mat')"