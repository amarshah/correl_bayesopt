function[ K12, K11invK12, K11invYminusmean, ...
          Kpred, Kpredinv, logdetKpred, mpred,...
          Sigma, Sigmainv, logdetSigma, lambda,...
          Kpredinvmpred, Sigmainvlambda, SigmaKpredinv ] = compute_pred_matrices( xstar, X, A, ls, Y, tau, mu, K11inv, mean, difs2_star )

[N,D] = size(X);
R     = size(A,1);

[ K12, K22 ] = compute_kernel_partial( X, xstar, A, ls, difs2_star );

Yminusmean        = Y - repmat(mean', N, 1);
K11invYminusmean  = K11inv*Yminusmean;
K11invK12         = K11inv*K12;
Kpred             = K22 - K12'*K11invK12;

C            = chol(Kpred + 1e-5*diag(diag(Kpred)));
Kpredinv     = solve_chol(C,eye(R));
logdetKpred  = 2*sum(log(diag(C)));

mpred     = K11invK12'*Yminusmean + mean';  % R x 1

Sigmainv     = Kpredinv + diag(1./tau);
C            = chol(Sigmainv + 1e-5*diag(diag(Sigmainv)));
Sigma        = solve_chol(C,eye(R));
logdetSigma  = -2*sum(log(diag(C)));

Kpredinvmpred = Kpredinv*mpred;
lambda        = Sigma*(Kpredinvmpred + mu./tau);

Sigmainvlambda = Kpredinvmpred + mu./tau;
SigmaKpredinv  = Sigma*Kpredinv;




