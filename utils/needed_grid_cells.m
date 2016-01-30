function[ C ] = needed_grid_cells( pareto,lower,upper )

Np    = size(pareto,1);
cells = all_grid_cells( [lower;pareto;upper] );

C = {};
for s=1:length(cells)
    % check if cell is inside pareto frontier
    corner = cells{s}(1,:);
    ind    = max(prod(double(bsxfun(@times,ones(Np,1),corner)<pareto),2));
    if ind==0
        C{end+1} = cells{s};
    end
end




