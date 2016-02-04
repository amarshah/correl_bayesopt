%%%%%%%%%%%
% A) load data from different files, 
% B) compute mean performance and standard devs 
% C) plot
%%%%%%%%%%%

task    = 'llvm';
methods = {'correl',  'multitask', 'indep', 'parego', 'random'};
colors  = { [0.9 0 0.1], ...      % red
            [0 0.8 0.2], ...      % green
            [1,0.8,0.3], ...      % yellow
            [0.5,0.6,1], ...      % blue
            [0.85,0.35,0.8] };    % purple
% methods = {'correl', 'indep', 'multitask', 'parego', 'random'};
% colors  = { [0.9 0 0.1], ...      % red
%             [0 0.8 0.2], ...      % green
%             [1,0.8,0.3], ...      % yellow
%             [0.5,0.6,1], ...      % blue
%             [0.85,0.35,0.8]};     % purple

%%%%%%   A   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = {};
for m=2:length(methods)
    pvs = [];
    for i=1:50
        method = methods{m};
        s = sprintf('./%s/%s/v%d.mat', task, method, i);        
        if exist(s, 'file')
            load(s);
            pvs(end+1,:) = PV;
        end        
    end 
    data{m} = pvs;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%   B   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
means = {};
medians = {};
stds = {};
for m=2:length(methods)
    means{m}   = mean(data{m});
    medians{m} = median(data{m});
    stds{m}    = std(data{m});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%   C   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
xlabel('t')
ylabel('Pareto Hypervolume')
pv = [];
for m=2:length(methods)
    pv(m) = plot(means{m}(1:100), 'LineWidth', 1.2, 'Color', colors{m});
    for i=1:100
       low  = means{m}(i) - stds{m}(i); 
       high = means{m}(i) + stds{m}(i); 
       if mod(i,4)==0
           plot([i i],[low high],'LineWidth',1,'Color', colors{m});
       end
    end    
end

r1=get(gca,'xtick');
r2=get(gca,'ytick');
gridxy(r1(2:end),r2(2:end),'color',[.8 .8 .8],'linewidth',0.5)

%legend(pv,'CEIPV-SLF','CEIPV-MT','IEIPV','ParEGO','Random','Location', 'SouthEast')
% matlab2tikz('.tikz','height','75mm','width','55mm')
% matlab2tikz('oka2.tikz','height','65mm','width','55mm')

 