#!/bin/sh
#
#$ -S /bin/bash
#
/usr/local/apps/matlab/matlabR2011b/bin/matlab - nodisplay -r "main( 'random', 'gp_threequarters.mat', 5, 8, 100, '/bigsratch/as793/multi_objective_bayesopt/gp_threequarters/random/v29.mat')