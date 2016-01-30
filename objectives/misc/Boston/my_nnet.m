% xtrain and xest are written in columns, so ytest and ytrain are row
% vectors

function[ error ] = my_nnet( xtrain, ytrain, xtest, ytest, W, b, weight_decay, penalty, n_iter )

N = size(xtrain,2);
L = length(W);
layer_sizes = zeros(L,1);
layer_sizes(1) = size(xtrain,1);
for l=2:L
    layer_sizes(l) = size(W{l},1);    
end

% % initialize Ws and bs
% W={}; b={};
% for l=2:h_layers+2
%     W{l} = randn(layer_sizes(l),layer_sizes(l-1))/sqrt(layer_sizes(l)*layer_sizes(l-1));
%     b{l} = zeros(layer_sizes(l),1);
% end


h = {xtrain};
dhdW = {}; dhdb = {};
for n=1:N
    dhdW{n} = {};
end

for iter=1:n_iter
    for l = 2:L
        if l==L
            h{l} = W{l}*h{l-1}+repmat(b{l},1,N);           % layer_sizes(l) x N
        else
            h{l} = max(W{l}*h{l-1}+repmat(b{l},1,N),0);    % layer_sizes(l) x N            
        end

        for n=1:N
            dhdhprev = bsxfun(@times,h{l}(:,n)>0,W{l});   % layer_sizes(l) x layer_sizes(l-1)   

            if l==L
                dhdW{n}{l} = reshape(h{l-1}(:,n),1,1,layer_sizes(l-1));
                dhdb{n}{l} = ones(layer_sizes(l-1),1);
            else
                tempW = zeros(layer_sizes(l),layer_sizes(l),layer_sizes(l-1));
                for i=1:layer_sizes(l)
                    for j=1:layer_sizes(l-1)
                        tempW(i,i,j) = (h{l}(i,n)>0)*h{l-1}(j,n);
                    end
                end
                dhdW{n}{l} = tempW;            
                dhdb{n}{l} = diag(h{l}(:,n));
            end

            if l>2
                for p=2:l-1                    
                    [~,I,J] = size(dhdW{n}{p});
                    tempW = zeros(layer_sizes(l),I,J);
                        for j=1:J
                            tempW(:,:,j) = reshape( dhdhprev*dhdW{n}{p}(:,:,j),layer_sizes(l),I,1 );
                        end
                    dhdW{n}{p} = tempW;
                    dhdb{n}{p} = dhdhprev*dhdb{n}{p};
                end
            end
        end        
    end

    % update each W and b
    for l=2:L
        dW = penalty*W{l};
        db = penalty*b{l};
        for n=1:N
            dW = dW - (ytrain(n)-h{L}(n))*reshape(dhdW{n}{l}(1,:,:),layer_sizes(l),layer_sizes(l-1)); 
            db = db - (ytrain(n)-h{L}(n))*dhdb{n}{l}(1,:)'; 
        end
        
        W{l} = W{l} - weight_decay*dW;
        b{l} = b{l} - weight_decay*db;
    end
end


Ntest = size(xtest,2);
pred  = xtest;
for l = 2:L-1
    pred = max(W{l}*pred+repmat(b{l},1,Ntest),0);
end
pred = W{L}*pred+repmat(b{L},1,Ntest);

error = sum((ytest-pred).*(ytest-pred))/Ntest;    

    
    
