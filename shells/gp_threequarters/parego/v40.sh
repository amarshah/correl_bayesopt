#!/bin/sh
#
#$ -S /bin/bash
#
/usr/local/apps/matlab/matlabR2011b/bin/matlab - nodisplay -r "main( 'parego', 'gp_threequarters.mat', 5, 8, 100, '/bigscratch/as793/multi_objective_bayesopt/gp_threequarters/parego/v40.mat')"