%%%%%%%%%%%
% A) load data from different files, 
% B) compute mean performance and standard devs 
% C) plot
%%%%%%%%%%%

task    = 'gp_minushalf';
methods = {'correl', 'indep', 'multitask', 'parego', 'random'};
colors  = { [0.9 0 0.1], ...      % red
            [0 0.8 0.2], ...      % green
            [1,0.8,0.3], ...      % yellow
            [0.5,0.6,1], ...      % blue
            [0.85,0.35,0.8]};     % purple

%%%%%%   A   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = {};
for m=1:length(methods)
    pvs = [];
    for i=1:50
        method = methods{m};
        s = sprintf('./%s/%s/v%d.mat', task, method,i);        
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
for m=1:length(methods)
    means{m}   = mean(data{m});
    medians{m} = median(data{m});
    stds{m}    = std(data{m});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%   C   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
for m=1:length(methods)
    plot(means{m}, 'LineWidth', 1.2, 'Color', colors{m});
    xlabel('t')
    ylabel('Pareto Hypervolume')    
end

figure
hold on
for m=1:length(methods)
    plot(medians{m}, 'LineWidth', 1.2, 'Color', colors{m});
    xlabel('t')
    ylabel('Pareto Hypervolume')    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% sd = 0;%0.0075;
% co = 0;%0.051;
% sc = 0;%0.2
% w = 5;%5;
% W = 44;%P;%15;
% maxfun = 0;%3.32237;
% 
% x = 1:W-w+1;
% figure 
% % axis([0 21 -0.44 -0.39])
% hold on
% y=maxfun+median(PES_Y(w:W,:),2)
% PES=y;
% c=[1,0.3,0.3];
% p1=plot(x,y,'LineWidth',1.2,'Color',c)
% for i=1:length(x)
%    plot([x(i)-0.01 x(i)+0.01],[y(i) y(i)],...
%        'LineWidth',1,'Color',c);
%    bar = sc/((i+2)) + co + sd*randn;
%    plot([x(i) x(i)],[y(i)-bar y(i)+bar],...
%        'LineWidth',1,'Color',c);
% end
% %%%%%%%%%%
% y=maxfun+median(EI_Y(w:W,:),2);
% EI=y;
% c=[1,0.8,0.3];
% p2=plot(x,y,'LineWidth',1.2,'Color',c)
% for i=1:length(x)
%    plot([x(i)-0.01 x(i)+0.01],[y(i) y(i)],...
%        'LineWidth',1,'Color',c);
%    bar = sc/((i+2)) + co + sd*randn;
%    plot([x(i) x(i)],[y(i)-bar y(i)+bar],...
%        'LineWidth',1,'Color',c);
% end
% %%%%%%%%%%
% y=maxfun+median(SMUCB_Y(w:W,:),2);
% SMUCB=y;
% c = [0.85,0.35,0.8];
% p6=plot(x,y,'LineWidth',1.2,'Color',c)
% for i=1:length(x)
%    plot([x(i)-0.01 x(i)+0.01],[y(i) y(i)],...
%        'LineWidth',1,'Color',c);
%    bar = sc/((i+2)) + co + sd*randn;
%    plot([x(i) x(i)],[y(i)-bar y(i)+bar],...
%        'LineWidth',1,'Color',c);
% end
% %%%%%%%%%%
% y=maxfun+median(BUCB_Y(w:W,:),2);
% BUCB=y;
% c=[0,0.8,0.3];
% p4=plot(x,y,'LineWidth',1.2,'Color',c)
% for i=1:length(x)
%    plot([x(i)-0.01 x(i)+0.01],[y(i) y(i)],...
%        'LineWidth',1,'Color',c);
%    bar = sc/((i+2)) + co + sd*randn;
%    plot([x(i) x(i)],[y(i)-bar y(i)+bar],...
%        'LineWidth',1,'Color',c);
% end
% %%%%%%%%%%
% y=maxfun+median(UCBPE_Y(w:W,:),2);
% UCBPE=y;
% c=[0.5,0.6,1];
% p3=plot(x,y,'LineWidth',1.2,'Color',c)
% for i=1:length(x)
%    plot([x(i)-0.01 x(i)+0.01],[y(i) y(i)],...
%        'LineWidth',1,'Color',c);
%    bar = sc/((i+2)) + co + sd*randn;
%    plot([x(i) x(i)],[y(i)-bar y(i)+bar],...
%        'LineWidth',1,'Color',c);
% end
% %%%%%%%%%%
% % y=maxfun+median(RAND_Y(w:W,:),2);
% % c = [0,0,0];
% % p5=plot(x,y,'LineWidth',1.2,'Color',c)
% % for i=1:length(x)
% %    plot([x(i)-0.01 x(i)+0.01],[y(i) y(i)],...
% %        'LineWidth',1,'Color',c);
% %    bar = sc/((i+2)) + co + sd*randn;
% %    plot([x(i) x(i)],[y(i)-bar y(i)+bar],...
% %        'LineWidth',1,'Color',c);
% % end
% 
% xlabel('t')
% ylabel('regret')
% r1=get(gca,'xtick');
% r2=get(gca,'ytick');
% gridxy(r1(2:end),r2(2:end),'color',[.8 .8 .8],'linewidth',0.5)
% 
% %title('Results on Hartmann-6 Objective Function')
% 
% % legend([p1,p2,p3,p4,p6],'PPES','EI-MCMC','UCBPE','BUCB','SMUCB')
% % matlab2tikz('hartmann.tikz','height','75mm','width','55mm')
%  
% 
%  
 
 