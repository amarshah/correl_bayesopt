% points are   NpxR pareto points and upper and lower points

function[ C ] = all_grid_cells( points )

[Np,R]   = size(points);
ordered  = sort(points);

C = {};
if R>1
    cells = all_grid_cells(points(:,2:R));
    for s=1:length(cells)
        for n=1:Np-1
            C{end+1} = [[ordered(n);ordered(n+1)],cells{s}];
        end    
    end
else
    for n=1:Np-1
        C{n} = [ordered(n);ordered(n+1)];
    end
end





