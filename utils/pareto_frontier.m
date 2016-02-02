% function takes outputs and returns a list of pareto optimal points
% Y is NxR 
function[ pareto_points ] = pareto_frontier( Y ) 

[N,R] = size(Y);

pareto_points = [];
for n=1:N
    y  = Y(n,:);
    
    Yt      = Y;
    Yt(n,:) = [];

    ind = bsxfun(@times,ones(N-1,1),y) <= Yt;
    ind = prod(double(ind),2);
    
    if sum(ind)==0
        pareto_points(end+1,:) = y;
    end
end

% remove repeats in pareto_points
Np = size(pareto_points, 1);
for np=1:Np-1
    p_point = pareto_points(np, :);
    difs2   = repmat(p_point,Np-np,1) - pareto_points(np+1:Np, :);
    difs2   = sum(difs2.*difs2, 2);
    if min(difs2) < 1e-6
        pareto_points(np, :) = [];
    end
end

