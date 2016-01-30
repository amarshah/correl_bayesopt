clear
% input space
Xall = 0:0.02:1;

% objectives
g1 = @(x) (x+0.01).*(x-0.1).*(x-0.33).*(x-0.65).*(x-0.95).*(x-1.5).*(x-0.73).*(x-1.1).*sin(x)*8500.*exp(-x);
f1 = @(x) (x+0.04).*(x-0.1).*(x-0.33).*(x-0.65).*(x-0.95).*(x-1.5).*cos(x*1.5)*560 + 0.1;
f2 = @(x) (x+0.03).*(x-0.17).*(x-0.4).*(x-0.72).*(x-1)*180+0.2;
f3 = @(x) f1(x+0.1).*sin(x-0.05).*(x-1.2)*3.3 + 0.12;
f4 = @(x) f3(x+0.3).*(x-1.1);
f5 = @(x) -0.7*g1(x) - 0.6.*f2(x);
obj = @(x) [g1(x-0.015), f5(x-0.015)]';

lower = [-1, -1];
upper = [ 1, 1];
xmin = [0];
xmax = [1];

nSamples = 20;
nIter = 1;

% observations
x = [0.07, 0.14, 0.41, 0.612, 0.825, 0.95, 1]';
y = obj(x)';

plot(Xall, g1(Xall-0.015), Xall, f5(Xall-0.015))
hold on 
scatter(x, y(:,1), 'x', 'LineWidth', 1.5, 'MarkerEdgeColor', 'k')
scatter(x, y(:,2), 'x', 'LineWidth', 1.5, 'MarkerEdgeColor', 'k')
%scatter(x, y(:,3), 'x', 'LineWidth', 1.5, 'MarkerEdgeColor', 'k')
hold off
init_y = y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute actual increase in vol %%%%%%%%%

IV = zeros(size(Xall));
F_cells = obs_to_integration_cells( init_y, lower, upper );    
cur_vol = compute_pareto_volume( F_cells, lower, upper );
for i=1:length(Xall)
    new = obj(Xall(i))';
    new_y = [init_y; new];
    F_cells = obs_to_integration_cells( new_y, lower, upper );    
    new_vol = compute_pareto_volume( F_cells, lower, upper );
    IV(i) = new_vol-cur_vol;
end

%%%%%% build curve for independent %%%%%%%%
D = length(xmin);
R = length(lower);

% sample random points 
N = size(x,1);
X = x;
Y = zeros(N,R);
F = zeros(N,R);
for nI=1:N
    y = obj(X(nI,:));
    Y(nI,:) = y';
    f = obj(X(nI,:));
    F(nI,:) = f';
end

% compute difX2
difX2 = zeros(N+1,N+1,D);
for n=1:N
    temp = bsxfun(@times,ones(N-n,1),X(n,:)) - X(n+1:N,:);     % N-n x D    
    difX2(n,n+1:N,:) = reshape(temp.*temp,1,N-n,D);
end 

% initialise samples
amp_samples = {};
ls_samples = {};
nvar_samples = {};
mean_samples = {};

v = Y-repmat(mean(Y),N,1);
v = sum(v.*v)/(N-1);      %1xR
 amp_samples{1} = 0.95*v;      %RxR
        temp    = 0.8*sqrt(xmax-xmin); %1xD
  ls_samples{1} = repmat(temp,R,1);  %RxD
nvar_samples{1} = 0.03*v;
mean_samples{1} = mean(Y);

%   sample hyperparams                      
[ amp_samples, ls_samples, nvar_samples, Kinv_samples, mean_samples, obj_means ] = sample_hypers_indep( amp_samples{end}, ls_samples{end}, nvar_samples{end}, mean_samples{end}, X, difX2, Y, nSamples, lower, upper );     

%   find integration cells based on prediction of noiseless objectives             
[ int_cells, ~ ] = obs_to_integration_cells( Y, lower, upper );

%   optimize EIV                
neg_acqui = @(x) neg_log_EIV_indep_all_cells_all_samples( x, X, Y, amp_samples, ls_samples, nvar_samples, Kinv_samples, mean_samples, int_cells );

log_EIV_indep = zeros(size(Xall));
for i=1:length(Xall)
    log_EIV_indep(i) = -neg_acqui(Xall(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% build curve for correl %%%%%%%%%%%

% initialise samples
v = Y-repmat(mean(Y),N,1);
v = sum(v.*v)/(N-1);      %1xR
   A_samples{1} = eye(R);      %RxR
        temp    = 0.8*sqrt(xmax-xmin); %1xD
  ls_samples{1} = repmat(temp,R,1);  %RxD
nvar_samples{1} = 0.02*v;
mean_samples{1} = mean(Y);

%   sample hyperparams                      
[ A_samples, ls_samples, nvar_samples, Kinv_samples, mean_samples, obj_means ] = sample_hypers_correl( A_samples{end}, ls_samples{end}, nvar_samples{end}, mean_samples{end}, X, difX2, Y, nSamples, lower, upper ); 

%   find integration cells based on prediction of noiseless objectives             
[ int_cells, ~ ] = obs_to_integration_cells( Y, lower, upper );

%   optimize EIV                
neg_acqui = @(x) neg_log_EIV_approx_all_cells_all_samples( x, X, Y, A_samples, ls_samples, nvar_samples, Kinv_samples, mean_samples, int_cells );

log_EIV_correl = zeros(size(Xall));
for i=1:length(Xall)
    log_EIV_correl(i) = -neg_acqui(Xall(i));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% log_EIV_correl_MCMC = zeros(size(Xall));
% for i=1:length(Xall)
%     log_EIV_correl_MCMC(i) = neg_log_EIV_MCMC_all_cells_all_samples( Xall(i), X, Y, A_samples, ls_samples, nvar_samples, Kinv_samples, mean_samples, int_cells );
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(Xall, IV, 'k')

figure
plot(Xall, log_EIV_indep, 'r', Xall, log_EIV_correl, 'g')
% figure
% plot(Xall, log_EIV_correl_MCMC, 'b')



