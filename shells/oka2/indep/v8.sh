#!/bin/sh
#
#$ -S /bin/bash
#
/usr/local/apps/matlab/matlabR2013a/bin/matlab - nodisplay -r "main( 'indep', 'oka2.mat', 5, 8, 100, '/bigscratch/as793/multi_objective_bayesopt/oka2/indep/v8.mat')"