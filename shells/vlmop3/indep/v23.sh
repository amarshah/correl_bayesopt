#!/bin/sh
#
#$ -S /bin/bash
#
/usr/local/apps/matlab/matlabR2013a/bin/matlab - nodisplay -r "main( 'indep', 'vlmop3.mat', 5, 8, 100, '/bigscratch/as793/multi_objective_bayesopt/vlmop3/indep/v23.mat')"