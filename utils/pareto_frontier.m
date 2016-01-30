% function takes outputs and returns a list of pareto optimal points
% Y is NxR 
function[ pareto_points ] = pareto_frontier( Y ) 

[N,R] = size(Y);

pareto_points = [];
for n=1:N
    y  = Y(n,:);
    
    Yt      = Y;
    Yt(n,:) = [];

    ind = bsxfun(@times,ones(N-1,1),y) < Yt;
    ind = prod(double(ind),2);
    
    if sum(ind)==0
        pareto_points(end+1,:) = y;
    end
end



