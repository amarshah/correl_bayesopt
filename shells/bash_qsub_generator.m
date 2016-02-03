
tasks  = {'boston'};
%{'oka2','dtlz1a','gp_minushalf','gp_threequarters'};
models  = {'indep','multitask'};
%{'indep','correl','multitask','parego','random'};


for t=1:length(tasks)
    task = tasks{t};
    bash_file = sprintf('../%s_bash.sh', task);
    
    fid = fopen(bash_file,'wt');
    for expt=1:50
        for m=1:length(models)
            model = models{m};
            unix = sprintf('./shells/%s/%s/%d.sh', task, model, expt);
        
            part1 = 'qsub ';
            part2 = sprintf('-o /bigscratch/as793/multi_objective_bayesopt/%s/%s/out%d ', task, model, expt);
            part3 = sprintf('-e /bigscratch/as793/multi_objective_bayesopt/%s/%s/error%d ', task, model, expt);
            part4 = sprintf('./shells/%s/%s/v%d.sh', task, model, expt);
            
            fprintf(fid,'%s',part1);
            fprintf(fid,'%s',part2);
            fprintf(fid,'%s',part3);
            fprintf(fid,'%s\n',part4);
        end
    end
    fclose(fid);
end