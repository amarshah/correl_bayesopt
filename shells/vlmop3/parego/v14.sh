#!/bin/sh
#
#$ -S /bin/bash
#
/usr/local/apps/matlab/matlabR2011b/bin/matlab - nodisplay -r "main( 'parego', 'vlmop3.mat', 5, 8, 100, '/bigscratch/as793/multi_objective_bayesopt/vlmop3/parego/v14.mat')"