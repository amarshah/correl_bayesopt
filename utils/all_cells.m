% returns all cells from a pareto set
% B is (Np+2)xR

function[ cells ] = all_cells( B )

[N,R] = size(B);

if R>0
    sub_cells = all_cells(B(:,2:end));

    b = sort(B(:,1));
    cells = {};
    S = length(sub_cells);
    for n=1:N-1
        if S>0
            for s=1:S
                temp  = [[b(n);b(n+1)],sub_cells{s}];
                cells = [cells,temp];
            end
        else
            cells = [cells,[b(n);b(n+1)]];
        end
    end
else
    cells = {};    
end


