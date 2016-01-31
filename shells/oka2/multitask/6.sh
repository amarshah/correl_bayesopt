#!/bin/sh
#
#$ -S /bin/bash
#
/usr/local/apps/matlab/matlabR2011b/bin/matlab - nodisplay -r "main( 'multitask', 'oka2.mat', 5, 8, 100, '/bigsratch/as793/multi_objective_bayesopt/oka2/multitask/6.mat')