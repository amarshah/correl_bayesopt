function[ y ] = boston_nnet( x )

load boston_nnet_params

weight_decay = exp( log(1e-5)+x(1)*(log(1e-2)-log(1e-5)) );
n_iter = round(10+40*x(2));
layer_size = 5 + round(20*x(3)); 

L1=13; L2=layer_size; L3=10; L4=1;

W={[]};
W{2} = (rand(L2,L1)-0.5)*2*sqrt(6/(L1+L2));
W{3} = (rand(L3,L2)-0.5)*2*sqrt(6/(L2+L3));
W{4} = (rand(L4,L3)-0.5)*2*sqrt(6/(L3+L4));

b={[]};
b{2} = zeros(L2,1);
b{3} = zeros(L3,1);
b{4} = [0];

y = my_nnet( Xtrain, Ytrain, Xtest, Ytest, W, b, weight_decay, penalty, n_iter );

y = [1 - min(y,2.5)/2.5, 1-x(3)*x(2)]';




