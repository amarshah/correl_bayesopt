#!/bin/sh
#
#$ -S /bin/bash
#
/usr/local/apps/matlab/matlabR2013a/bin/matlab - nodisplay -r "main( 'multitask', 'llvm.mat', 5, 8, 100, '/bigscratch/as793/multi_objective_bayesopt/llvm/multitask/v31.mat')"