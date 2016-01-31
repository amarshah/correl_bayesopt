#!/bin/sh
#
#$ -S /bin/bash
#
/usr/local/apps/matlab/matlabR2011b/bin/matlab - nodisplay -r "main( 'parego', 'dtlz1a.mat', 5, 8, 100, '/bigsratch/as793/multi_objective_bayesopt/dtlz1a/parego/v11.mat')