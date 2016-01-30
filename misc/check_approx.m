%%% check that approx gives similar acquisition value to MCMC %%%
clear
rng shuffle

f1 = @(x) cos(6*x(1)*x(2)-0.5)*sqrt(0.1+x(1)) + x(2)/4;
f2 = @(x) sin(4*(0.75-x(1))*(x(2)+0.25))*sqrt(0.8+x(2)) + sqrt(x(1))-x(1);
obj = @(x) [f1(x), f2(x)]';

obj1 = zeros(ng,ng);
obj2 = zeros(ng,ng);
for i=1:ng
    for j=1:ng
        x = [Xgrid(i), Xgrid(j)];
        obj1(i,j) = f1(x);
        obj2(i,j) = f2(x);
    end
end


xmin = [0,0];
xmax = [1,1];

lower = [-1.5, -1.5];
upper = [ 2, 2];

nSamples = 10;
nInitial = 10;

Xinit = rand(nInitial, length(xmin));
Xgrid = 0:0.02:1;
ng = length(Xgrid);

%%%%%%%%%%%%%%

D = length(xmin);
R = length(lower);
N = nInitial;
X = Xinit;
Y = zeros(N,R);
for nI=1:N
    y = obj(X(nI,:));
    Y(nI,:) = y';
end

% compute difX2
difX2 = zeros(N+1,N+1,D);
for n=1:N
    temp = bsxfun(@times,ones(N-n,1),X(n,:)) - X(n+1:N,:);     % N-n x D    
    difX2(n,n+1:N,:) = reshape(temp.*temp,1,N-n,D);
end 

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
acqui = @(x) - neg_log_EIV_approx_all_cells_all_samples( x, X, Y, A_samples, ls_samples, nvar_samples, Kinv_samples, mean_samples, int_cells );
acqui_MCMC = @(x) - neg_log_EIV_MCMC_all_cells_all_samples( x, X, Y, A_samples, ls_samples, nvar_samples, Kinv_samples, mean_samples, int_cells );

log_EIV_correl = zeros(ng, ng);
log_EIV_correl_MCMC = zeros(ng, ng);
for i=1:ng
    for j=1:ng
        x = [Xgrid(i), Xgrid(j)];
        log_EIV_correl(i, j) = acqui(x);
    end
end

for i=1:ng
    i
    for j=1:ng
        x = [Xgrid(i), Xgrid(j)];
        log_EIV_correl_MCMC(i, j) = acqui_MCMC(x);
    end
end


%%%%%%%%%%%%%%%%%%%%%


mymap=[255,255,255;
253,255,255;
251,255,255;
249,255,255;
247,255,255;
245,255,255;
243,255,255;
241,255,255;
239,255,255;
236,255,255;
234,255,255;
232,255,255;
230,255,255;
228,255,255;
226,255,255;
224,255,255;
222,255,255;
220,255,255;
218,255,255;
216,255,255;
214,255,255;
212,255,255;
210,255,255;
208,255,255;
206,255,255;
203,255,255;
201,255,255;
199,255,255;
197,255,255;
195,255,255;
193,255,255;
191,255,255;
189,255,255;
187,255,255;
185,255,255;
183,255,255;
181,255,255;
179,255,255;
177,255,255;
175,255,255;
173,255,255;
171,255,255;
168,255,255;
166,255,255;
164,255,255;
162,255,255;
160,255,255;
158,255,255;
156,255,255;
154,255,255;
152,255,255;
150,255,255;
148,255,255;
146,255,255;
144,255,255;
142,255,255;
140,255,255;
138,255,255;
135,255,255;
133,255,255;
131,255,255;
129,255,255;
127,255,255;
125,255,255;
123,254,255;
121,252,255;
119,250,255;
117,248,255;
115,246,255;
113,244,255;
112,242,255;
110,240,255;
108,238,255;
106,236,255;
104,234,255;
102,232,255;
100,230,255;
98,228,255;
96,226,255;
94,224,255;
92,222,255;
90,219,255;
88,217,255;
86,215,255;
84,213,255;
82,211,255;
80,209,255;
79,207,255;
77,205,255;
75,203,255;
73,201,255;
71,199,255;
69,197,255;
67,195,255;
65,193,255;
63,191,255;
61,189,255;
59,186,255;
57,184,255;
55,182,255;
53,180,255;
51,178,255;
49,176,255;
48,174,255;
46,172,255;
44,170,255;
42,168,255;
40,166,255;
38,164,255;
36,162,255;
34,160,255;
32,158,255;
30,156,255;
28,154,255;
26,151,255;
24,149,255;
22,147,255;
20,145,255;
18,143,255;
16,141,255;
15,139,255;
13,137,255;
11,135,255;
9,133,255;
7,131,255;
5,129,255;
3,127,255;
1,125,255;
0,123,254;
0,121,252;
0,119,250;
0,117,248;
0,115,246;
0,113,244;
0,111,242;
0,109,240;
0,107,237;
0,105,235;
0,103,233;
0,101,231;
0,99,229;
0,97,227;
0,96,225;
0,94,223;
0,92,221;
0,90,219;
0,88,217;
0,86,215;
0,84,213;
0,82,211;
0,80,209;
0,78,207;
0,76,205;
0,74,202;
0,72,200];
 
mymap2 = mymap;
mymap2(:,3) = mymap(:,2);
mymap2(:,2) = mymap(:,3);

figure
hold on
contourf(linspace(0,1,51),...
         linspace(0,1,51),...
         0.4*exp(log(3)+temp)+0.6*exp(log_EIV_correl_MCMC),...
         8,'LineStyle','none')
colormap(mymap/255)
cb=colorbar;
scatter(Xinit(:,1),Xinit(:,2),'k','s','LineWidth',0.3,'MarkerFaceColor','k')
plot([0 1],[0 0],'k')
plot([0 0],[0 1],'k')
plot([0 1],[1 1],'k')
plot([1 1],[0 1],'k')
hold off
matlab2tikz('EIV_approx.tikz','height','75mm','width','75mm')

figure
hold on
contourf(linspace(0,1,51),...
         linspace(0,1,51),...
         exp(log_EIV_correl_MCMC),...
         8,'LineStyle','none')
colormap(mymap/255)
cb=colorbar;
scatter(Xinit(:,1),Xinit(:,2),'k','s','LineWidth',0.3,'MarkerFaceColor','k')
plot([0 1],[0 0],'k')
plot([0 0],[0 1],'k')
plot([0 1],[1 1],'k')
plot([1 1],[0 1],'k')
hold off
matlab2tikz('EIV_mcmc.tikz','height','75mm','width','75mm')


figure
hold on
contourf(linspace(0,1,51),...
         linspace(0,1,51),...
         obj1*1.3,...
         8,'LineStyle','none')
colormap(mymap2/255)
cb=colorbar;
scatter(Xinit(:,1),Xinit(:,2),'k','s','LineWidth',0.3,'MarkerFaceColor','k')
plot([0 1],[0 0],'k')
plot([0 0],[0 1],'k')
plot([0 1],[1 1],'k')
plot([1 1],[0 1],'k')
hold off
matlab2tikz('obj1.tikz','height','37.5mm','width','75mm')

figure
hold on
contourf(linspace(0,1,51),...
         linspace(0,1,51),...
         obj2,...
         8,'LineStyle','none')
colormap(mymap2/255)
cb=colorbar;
scatter(Xinit(:,1),Xinit(:,2),'k','s','LineWidth',0.3,'MarkerFaceColor','k')
plot([0 1],[0 0],'k')
plot([0 0],[0 1],'k')
plot([0 1],[1 1],'k')
plot([1 1],[0 1],'k')
hold off
matlab2tikz('obj2.tikz','height','37.5mm','width','75mm')


%save('check_approx.mat','obj1','obj2','Xinit','Xgrid','EIV_correl_approx','EIV_correl_numint')

