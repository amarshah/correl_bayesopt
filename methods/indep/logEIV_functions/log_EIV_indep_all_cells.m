% Y is NxR

function[g,dg] = log_EIV_indep_all_cells( xstar, X, Y, amp, ls, nvar, K22invs, mean, cells )

[N,D]  = size(X);
 R     = size(Y,2);
ncells = length(cells);

f  = 0;
df = zeros(length(xstar),1);

 m_preds = zeros(R*ncells,1); % R batches of ncellsx1
 v_preds = zeros(R*ncells,1);
dm_preds = zeros(R*ncells,D);
dv_preds = zeros(R*ncells,D);

da = zeros(R*ncells,D);
db = zeros(R*ncells,D);

cdf_a    = zeros(R*ncells,1);
cdf_b    = zeros(R*ncells,1);
pdf_a    = zeros(R*ncells,1);
pdf_b    = zeros(R*ncells,1);

a      = zeros(R*ncells,1);
b      = zeros(R*ncells,1);
lrs    = zeros(R*ncells,1);
urs    = zeros(R*ncells,1);

% R,R,DxR,DxR
[ m_pred_all, v_pred_all, dm_pred_all, dv_pred_all ] = compute_pred_matrices_small( xstar, X, Y, amp, ls, K22invs, mean );

for r=1:R
%     Yr = Y(:,r);
%     [ m_pred, v_pred, dm_pred, dv_pred ] = compute_pred_matrix_small( xstar, X, Yr, amp(r), ls(r,:), K22inv{r} );
    m_pred =  m_pred_all(r);
    v_pred =  v_pred_all(r);
   dm_pred = dm_pred_all(:,r);
   dv_pred = dv_pred_all(:,r);

     m_preds((r-1)*ncells+1:r*ncells) = repmat(m_pred,ncells,1);
     v_preds((r-1)*ncells+1:r*ncells) = repmat(v_pred,ncells,1);
    dm_preds((r-1)*ncells+1:r*ncells,:) = repmat(dm_pred',ncells,1);
    dv_preds((r-1)*ncells+1:r*ncells,:) = repmat(dv_pred',ncells,1);
    
    dm_pred_over_sqrt_v_pred        = dm_pred'/sqrt(v_pred);    %1xD
    dv_pred_over_v_pred             = dv_pred'/v_pred;          %1xD
    dm_pred_over_sqrt_v_pred_repeat = bsxfun(@times,ones(ncells,1),dm_pred_over_sqrt_v_pred);   % ncells x D
    
    for s=1:ncells
        lr = cells{s}(1,r); % lower(r);
        ur = cells{s}(2,r); % upper(r);
 
        lrs((r-1)*ncells+s) = lr;
        urs((r-1)*ncells+s) = ur;         
    end
     a((r-1)*ncells+1:r*ncells) = (urs((r-1)*ncells+1:r*ncells)-m_pred)/sqrt(v_pred);
     b((r-1)*ncells+1:r*ncells) = (lrs((r-1)*ncells+1:r*ncells)-m_pred)/sqrt(v_pred);
    
    da((r-1)*ncells+1:r*ncells,:) = -dm_pred_over_sqrt_v_pred_repeat - 0.5*bsxfun(@times,a((r-1)*ncells+1:r*ncells),dv_pred_over_v_pred);
    db((r-1)*ncells+1:r*ncells,:) = -dm_pred_over_sqrt_v_pred_repeat - 0.5*bsxfun(@times,b((r-1)*ncells+1:r*ncells),dv_pred_over_v_pred);

end

cdf_a = normcdf(a);
cdf_b = normcdf(b);
pdf_a = normpdf(a);
pdf_b = normpdf(b);

 temp = sqrt(v_preds).*(pdf_b-pdf_a) + (m_preds-lrs).*(cdf_a-cdf_b);        % R*ncells x 1
 temp = max(temp,1e-300);
%    f = sum(log(temp));     
other = m_preds-lrs+sqrt(v_preds);    
dtemp = -0.5*bsxfun(@times,temp./v_preds,dv_preds) + bsxfun(@times,cdf_a-cdf_b,dm_preds) ...
            - bsxfun(@times,other.*pdf_b,db) + bsxfun(@times,other.*pdf_a,da);    % R*ncells x D

%   df = sum(dtemp./repmat(temp,1,D))';     % D x 1
        
%%%%%%%%%%%%%%
 logEIV_cells = zeros(ncells, 1);
dlogEIV_cells = zeros(ncells, D);

for s=1:ncells   
     temp_cell =  temp(s:ncells:R*ncells);    % Rx1
    dtemp_cell = dtemp(s:ncells:R*ncells,:);  % RxD
    
     logEIV_cells(s)   = sum(log(temp_cell));  
    dlogEIV_cells(s,:) = sum(dtemp_cell./bsxfun(@times,temp_cell,ones(1,D)));  
end

 g = logsumexp(logEIV_cells);                % scalar
 p = exp(logEIV_cells - g);                  % ncells x 1
dg = sum(bsxfun(@times, p, dlogEIV_cells));  % 1 x D  



