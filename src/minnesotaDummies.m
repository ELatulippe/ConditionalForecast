function [Yd, Xd, info] = minnesotaDummies(Y, p, tcode, qExo, opts)
% MINNESOTADUMMIES  Minnesota / Normal-inverse-Wishart prior as dummy data.
%
%   [Yd, Xd] = minnesotaDummies(Y, p, tcode, qExo)
%   [Yd, Xd] = minnesotaDummies(Y, p, tcode, qExo, struct('lambda',0.1))
%
% Litterman's prior, implemented the way Banbura, Giannone and Reichlin (2010)
% do it: as artificial observations appended to the data. The posterior mean
% is then just OLS on the augmented sample, so nothing about the estimator
% changes -- only the design matrix gets taller. That is what keeps this a
% contained addition rather than a rewrite.
%
% The prior shrinks each equation toward a univariate random walk (for
% variables entering in levels) or toward white noise (for variables already
% differenced), with coefficients on longer lags shrunk harder. lambda sets
% the overall tightness: lambda -> Inf is OLS, lambda -> 0 pins the model at
% the prior. Values around 0.1-0.2 are typical for macro VARs of this size.
%
% BLOCKS (with the intercept FIRST in X, matching estimateVAR)
%   1. own and cross lags, scaled by lag number j and by sigma_i / lambda
%   2. a block that ties down Sigma_u
%   3. a near-flat prior on the intercept
% Exogenous regressors get all-zero columns, i.e. a flat prior -- pandemic and
% seasonal dummies should not be shrunk toward anything.
%
% INPUTS
%   Y      T x K data (the VAR's input, already transformed)
%   p      lag order
%   tcode  1 x K cellstr of 'dlog'/'level'; sets the prior mean on the first
%          own lag -- 1 for a level variable, 0 for one already differenced
%   qExo   number of exogenous regressors (0 if none)
%   opts   .lambda   overall tightness (default 0.2)
%          .epsilon  intercept looseness (default 1e-4; smaller = flatter)
%          .sigma    K x 1 prior scales; default = residual sd of a univariate
%                    AR(p) on each series
%          .delta    K x 1 prior mean on the first own lag; default from tcode
%
% OUTPUT
%   Yd, Xd  dummy blocks to append to the estimation sample
%   info    .lambda, .sigma, .delta, .nDummy

    [T, K] = size(Y);
    if nargin < 4 || isempty(qExo), qExo = 0; end
    if nargin < 5, opts = struct(); end
    if ~isfield(opts,'lambda')  || isempty(opts.lambda),  opts.lambda  = 0.2;  end
    if ~isfield(opts,'epsilon') || isempty(opts.epsilon), opts.epsilon = 1e-4; end

    % ---- prior scales: univariate AR(p) residual sd ----------------------
    if isfield(opts,'sigma') && ~isempty(opts.sigma)
        sigma = opts.sigma(:);
    else
        sigma = zeros(K,1);
        for i = 1:K
            y  = Y(p+1:T, i);
            Xi = ones(T-p, 1);
            for j = 1:p, Xi = [Xi, Y(p+1-j:T-j, i)]; end %#ok<AGROW>
            b  = Xi \ y;
            e  = y - Xi*b;
            sigma(i) = sqrt(e.'*e / (numel(e) - size(Xi,2)));
        end
    end
    if any(sigma <= 0) || any(~isfinite(sigma))
        error('minnesotaDummies:sigma', ...
            'Non-positive prior scale for variable(s) %s.', ...
            mat2str(find(sigma <= 0 | ~isfinite(sigma))'));
    end

    % ---- prior mean on the first own lag ---------------------------------
    if isfield(opts,'delta') && ~isempty(opts.delta)
        delta = opts.delta(:);
    else
        delta = zeros(K,1);
        for i = 1:K
            if strcmp(tcode{i}, 'level'), delta(i) = 1; end
        end
    end

    lam = opts.lambda;  eps0 = opts.epsilon;
    nX  = 1 + K*p + qExo;                  % [const, lags, exo]

    % ---- block 1: lag coefficients ---------------------------------------
    Yd1 = [diag(delta .* sigma) / lam; zeros(K*(p-1), K)];
    Xd1 = zeros(K*p, nX);
    Xd1(:, 2:1+K*p) = kron(diag(1:p), diag(sigma) / lam);

    % ---- block 2: scale of Sigma_u ---------------------------------------
    Yd2 = diag(sigma);
    Xd2 = zeros(K, nX);

    % ---- block 3: intercept ----------------------------------------------
    Yd3 = zeros(1, K);
    Xd3 = zeros(1, nX);  Xd3(1) = eps0;

    Yd = [Yd1; Yd2; Yd3];
    Xd = [Xd1; Xd2; Xd3];

    info = struct('lambda',lam, 'epsilon',eps0, 'sigma',sigma, ...
                  'delta',delta, 'nDummy',size(Yd,1));
end
