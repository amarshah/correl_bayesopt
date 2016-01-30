function[y] = boston_estimate(x)

load boston_values

y=0;

dif = repmat(x,size(Xgrid,1),1)-Xgrid;   %NxD
for i=1:5
    % compute K21
    temp  = 1./(lss(i,:)'.*lss(i,:)');   % Dx1
    K12   = amp(i)*exp(-0.5*(dif.*dif)*temp);   %Nx1 
    
    % predict
    pred = const_mean(i) + K12'*Kinvs{i}*(Ygrid-const_mean(i));
    
    y = y+pred/5;
end







