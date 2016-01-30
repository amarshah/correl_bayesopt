% Y is NxR

function[g] = log_EIV_indep_all_cells_noderiv( xstar, X, Y, amp, ls, nvar, K22invs, mean, cells )

[N,D]  = size(X);
 R     = size(Y,2);
ncells = length(cells);

f  = 0;

 m_preds = zeros(R*ncells,1); % R batches of ncellsx1
 v_preds = zeros(R*ncells,1);

cdf_a    = zeros(R*ncells,1);
cdf_b    = zeros(R*ncells,1);
pdf_a    = zeros(R*ncells,1);
pdf_b    = zeros(R*ncells,1);

a      = zeros(R*ncells,1);
b      = zeros(R*ncells,1);
lrs    = zeros(R*ncells,1);
urs    = zeros(R*ncells,1);

% R,R,DxR,DxR
[ m_pred_all, v_pred_all ] = compute_pred_matrices_small_noderiv( xstar, X, Y, amp, ls, K22invs, mean );

for r=1:R
    m_pred =  m_pred_all(r);
    v_pred =  v_pred_all(r);

     m_preds((r-1)*ncells+1:r*ncells) = repmat(m_pred,ncells,1);
     v_preds((r-1)*ncells+1:r*ncells) = repmat(v_pred,ncells,1);
        
    for s=1:ncells
        lr = cells{s}(1,r); % lower(r);
        ur = cells{s}(2,r); % upper(r);
 
        lrs((r-1)*ncells+s) = lr;
        urs((r-1)*ncells+s) = ur;         
    end
     a((r-1)*ncells+1:r*ncells) = (urs((r-1)*ncells+1:r*ncells)-m_pred)/sqrt(v_pred);
     b((r-1)*ncells+1:r*ncells) = (lrs((r-1)*ncells+1:r*ncells)-m_pred)/sqrt(v_pred);
    
end

cdf_a = normcdf(a);
cdf_b = normcdf(b);
pdf_a = normpdf(a);
pdf_b = normpdf(b);

 temp = sqrt(v_preds).*(pdf_b-pdf_a) + (m_preds-lrs).*(cdf_a-cdf_b);        % R*ncells x 1
 temp = max(temp,1e-300);

 logEIV_cells = zeros(ncells, 1);

for s=1:ncells   
     temp_cell =  temp(s:ncells:R*ncells);    % Rx1    
     logEIV_cells(s)   = sum(log(temp_cell));  
end

 g = logsumexp(logEIV_cells);                % scalar



