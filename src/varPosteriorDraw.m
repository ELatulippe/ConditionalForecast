function md = varPosteriorDraw(model)
% VARPOSTERIORDRAW  One draw from the VAR's normal-inverse-Wishart posterior.
%
%   md = varPosteriorDraw(model)      % model from estimateVAR(..., prior)
%
% With a conjugate prior implemented as dummy observations, the posterior is
% normal-inverse-Wishart in closed form on the augmented sample:
%
%     Sigma | Y      ~ IW( S, Tstar - nreg )
%     B | Sigma, Y   ~ MN( Bhat, Sigma kron inv(X*'X*) )
%
% so uncertainty is sampled directly rather than bootstrapped. That matters
% once a prior is in play: a residual bootstrap would resample around the
% shrunk point estimate while ignoring the prior's own contribution to
% posterior uncertainty, understating it in exactly the region -- long lags,
% weakly identified coefficients -- where the prior is doing the work.
%
% Returns a model struct of the same shape as estimateVAR's, so anything
% downstream (varMA, baselineForecast, irfCholesky) accepts it unchanged.
%
% Needs no toolbox: the inverse-Wishart draw is built from a Bartlett-style
% sum of outer products of normal vectors.

    if ~isfield(model,'S') || isempty(model.S)
        error('varPosteriorDraw:noPosterior', ...
            ['This model carries no posterior ingredients. Estimate it with ' ...
             'a prior: estimateVAR(Y, p, true, Xexo, prior).']);
    end

    K = model.K;  p = model.p;
    nu = model.Tstar - model.nreg;
    if nu <= K + 1
        error('varPosteriorDraw:df', ...
            'Only %d posterior degrees of freedom for K = %d.', nu, K);
    end

    % ---- Sigma ~ IW(S, nu) -----------------------------------------------
    % Draw W ~ Wishart(inv(S), nu) as sum of outer products, then invert.
    Sinv = inv(0.5*(model.S + model.S.'));
    Sinv = 0.5*(Sinv + Sinv.');
    Ls   = cholPSD(Sinv);
    Z    = Ls * randn(K, nu);
    W    = Z * Z.';
    Sig  = inv(0.5*(W + W.'));
    Sig  = 0.5*(Sig + Sig.');

    % ---- B | Sigma ~ MN(Bhat, Sigma kron inv(X'X)) -----------------------
    % vec(B - Bhat) = (P_Sig kron P_x) vec(Z2), so B = Bhat + P_x Z2 P_Sig'.
    Px   = cholPSD(0.5*(model.XtXinv + model.XtXinv.'));
    PSig = cholPSD(Sig);
    Z2   = randn(model.nreg, K);
    Beta = model.Beta + Px * Z2 * PSig.';

    % ---- repack in estimateVAR's shape -----------------------------------
    md = model;
    md.Beta   = Beta;
    md.Sigmau = Sig;

    Bt = Beta.';
    if model.includeConst
        md.nu = Bt(:,1);   Aall = Bt(:,2:end);
    else
        md.nu = zeros(K,1); Aall = Bt;
    end
    if model.qExo > 0
        md.Bexo = Aall(:, K*p+1 : K*p+model.qExo).';
        Aall    = Aall(:, 1:K*p);
    end
    for i = 1:p
        md.Acell{i} = Aall(:, (i-1)*K+1 : i*K);
    end
    A = zeros(K*p, K*p);
    A(1:K, :) = Aall;
    if p > 1, A(K+1:end, 1:K*(p-1)) = eye(K*(p-1)); end
    md.A  = A;
    md.mu = md.J.' * md.nu;
end

% ======================================================================
function L = cholPSD(M)
% Cholesky with a jitter fallback, for matrices that are positive definite in
% exact arithmetic but marginal in floating point.
    [L, flag] = chol(M, 'lower');
    if flag == 0, return; end
    d = mean(abs(diag(M)));
    for k = -12:2:0
        [L, flag] = chol(M + 10^k * d * eye(size(M,1)), 'lower');
        if flag == 0, return; end
    end
    error('varPosteriorDraw:chol', 'Matrix is not positive definite.');
end
