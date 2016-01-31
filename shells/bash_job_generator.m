
tasks  = {'oka2','dtlz1a','gp_threequarters','gp_minushalf'};
models = {'correl','indep','multitask','parego','random'};

line1  = '#!/bin/sh';
line2  = '#';
line3  = '#$ -S /bin/bash';
line4  = '#';
unixstr1 = '/usr/local/apps/matlab/matlabR2011b/bin/matlab - nodisplay -r ';

nInitial = 5;
nSamples = 8;
nIter    = 100;

for t=1:length(tasks)
    task  = tasks{t};
    for m=1:length(models)
        model = models{m};
        for expt=1:50
            
            savefile = sprintf('/bigsratch/as793/multi_objective_bayesopt/%s/%s/v%d.mat', task, model, expt);
            bash_save_file = sprintf('./%s/%s/v%d.sh', task, model, expt);
            unixstr2 = sprintf('"main( ''%s'', ''%s.mat'', %d, %d, %d, ''%s'')', model, task, nInitial, nSamples, nIter, savefile);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            fid = fopen(bash_save_file,'wt');

            fprintf(fid,'%s\n',line1);
            fprintf(fid,'%s\n',line2);
            fprintf(fid,'%s\n',line3);
            fprintf(fid,'%s\n',line4);
            fprintf(fid,'%s',unixstr1);
            fprintf(fid,'%s',unixstr2);

            fclose(fid);
        end
    end
end
