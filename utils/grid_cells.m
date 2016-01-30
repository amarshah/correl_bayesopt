% function takes pareto points and defines cells C for integration
% pareto   is Np x R

function[ C ] = grid_cells( pareto, lower, upper )

[Np,R] = size(pareto);

if Np>0
    B = [lower;pareto;upper];   
    sub_cells = all_cells( B(:,2:end) );
    [b,ind] = max(pareto(:,1));

    cells = {};
    for s=1:length(sub_cells)
        temp  = [[b;upper(1)],sub_cells{s}];
        cells = [cells,temp];
    end

    new_pareto        = pareto;
    new_pareto(ind,:) = [];
    new_upper         = upper;
    new_upper(1)      = b;
    if Np>1
        new_lower = min(pareto);
    else
        new_lower = pareto;
    end
    new_lower(1)      = lower(1);
    C = [cells,grid_cells(new_pareto,new_lower,new_upper)];
else
    C = {[lower;upper]};    
end

