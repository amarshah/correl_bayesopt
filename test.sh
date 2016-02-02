#!/bin/sh
#
#$ -S /bin/bash
#
/usr/local/apps/matlab/matlabR2011b/bin/matlab -nodisplay -r "main( 'indep', 'oka2.mat', 8, 8, 1, '/bigscratch/as793/test.mat')"
