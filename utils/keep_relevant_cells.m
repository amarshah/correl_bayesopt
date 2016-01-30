
function[ cells ] = keep_relevant_cells( C, pareto )

[Np,R] = size(pareto);

cells = {};
for s=1:length(C)
    upper = C{s}(2,:);
    
    if Np>1
        ind = bsxfun(@times,ones(Np,1),upper) - pareto;  %Np x R
    else
        ind = upper - pareto;
    end
    ind = max(double(ind'))';

    if min(ind)>0
        cells{end+1} = C{s};
    end
end



