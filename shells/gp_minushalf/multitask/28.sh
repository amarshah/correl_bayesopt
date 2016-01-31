#!/bin/sh
#
#$ -S /bin/bash
#
/usr/local/apps/matlab/matlabR2011b/bin/matlab - nodisplay -r "main( 'multitask', 'gp_minushalf.mat', 5, 8, 100, '/bigsratch/as793/multi_objective_bayesopt/gp_minushalf/multitask/28.mat')