#!/bin/sh
#
#$ -S /bin/bash
#
/usr/local/apps/matlab/matlabR2011b/bin/matlab - nodisplay -r "main( 'correl', 'dtlz1a.mat', 5, 8, 100, '/bigsratch/as793/multi_objective_bayesopt/dtlz1a/correl/v21.mat')