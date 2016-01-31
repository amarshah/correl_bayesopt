#!/bin/sh
#
#$ -S /bin/bash
#
/usr/local/apps/matlab/matlabR2011b/bin/matlab - nodisplay -r "main( 'parego', 'gp_minushalf.mat', 5, 8, 100, '/bigsratch/as793/multi_objective_bayesopt/gp_minushalf/parego/v12.mat')