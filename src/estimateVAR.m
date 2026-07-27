function model = estimateVAR(Y, p, includeConst, Xexo, prior)
% ESTIMATEVAR  OLS estimation of a reduced-form VAR(p) and its companion form.
%
%   model = estimateVAR(Y, p, includeConst)
%   model = estimateVAR(Y, p, includeConst, Xexo)
%   model = estimateVAR(Y, p, includeConst, Xexo, prior)   % Minnesota
%
% Implements Eqs. (1)-(2) of the derivation note.
%
% INPUTS
%   Y            T x K matrix of (already transformed / stationary) data,
%                rows = time, columns = variables, ordered as they enter y_t.
%   p            lag order.
%   includeConst (optional, default true) include intercept nu.
%   Xexo         (optional) T x q exogenous regressors aligned row-for-row
%                with Y, e.g. pandemic dummies from covidDummies. They enter
%                every equation, are estimated alongside nu and A_i, and are
%                appended AFTER the lag block so nu and Acell unpack
%                unchanged. Because such dummies are zero outside their
%                window, they are zero over the forecast horizon and the
%                baseline recursion is unaffected -- their whole role is to
%                keep those observations out of Sigma_u.
%   prior        (optional) struct with .Yd and .Xd from minnesotaDummies.
%                The dummy rows are appended to the estimation sample, so the
%                posterior mean is OLS on the augmented data -- the estimator
%                itself is unchanged. The posterior ingredients are returned
%                so varPosteriorDraw can sample from it:
%                  Sigma | Y ~ IW(S, Tstar - nreg)
%                  B | Sigma, Y ~ MN(Bhat, Sigma kron inv(X*'X*))
%
% OUTPUT: struct 'model' with fields
%   .nu      K x 1        intercept
%   .Acell   1 x p cell   {A_1,...,A_p}, each K x K
%   .A       Kp x Kp      companion matrix A          [Eq. (2)]
%   .mu      Kp x 1       companion intercept  mu = J'*nu
%   .J       K x Kp       selection matrix [I_K 0 ... 0]
%   .Sigmau  K x K        residual covariance  Sigma_u = U'U/(Teff - Kp - 1)
%   .U       Teff x K     OLS residuals (reduced-form innovations)
%   .Beta    (1+Kp+q) x K stacked OLS coefficients
%   .Bexo    q x K        exogenous-regressor coefficients (empty if none)
%   .qExo    scalar       number of exogenous regressors
%   .prior   struct        the prior info, empty when none was used
%   .XtXinv  k x k         inv(X*'X*) on the augmented sample
%   .S       K x K         posterior scale matrix (residual cross-product)
%   .Tstar   scalar        rows in the augmented sample
%   .nreg    scalar        regressors per equation
%   .K,.p,.T,.Teff        dimensions
%
% The variance normaliser follows the note: Teff - (K*p + 1).

    if nargin < 3 || isempty(includeConst), includeConst = true; end
    if nargin < 4, Xexo = []; end
    if nargin < 5, prior = []; end
    [T, K] = size(Y);
    if ~isempty(Xexo) && size(Xexo,1) ~= T
        error('estimateVAR:exoRows', ...
            'Xexo has %d rows but Y has %d; they must align.', size(Xexo,1), T);
    end
    Teff = T - p;

    % Build regressor matrix X (Teff x (1+Kp)) and response Yr (Teff x K).
    Yr = Y(p+1:T, :);
    X  = zeros(Teff, K*p);
    for i = 1:p
        X(:, (i-1)*K+1 : i*K) = Y(p+1-i : T-i, :);   % lag i
    end
    if includeConst
        X = [ones(Teff,1), X];
    end
    qExo = size(Xexo, 2);
    if qExo > 0
        X = [X, Xexo(p+1:T, :)];        % appended last: nu and A_i unpack as before
    end

    % Append the prior's dummy observations, if any. Everything below then
    % runs on the augmented sample, which is exactly what makes a Minnesota
    % prior a change of design matrix rather than a change of estimator.
    Tobs = size(X,1);
    if ~isempty(prior) && isfield(prior,'Yd') && ~isempty(prior.Yd)
        if size(prior.Xd,2) ~= size(X,2)
            error('estimateVAR:priorCols', ...
                'Prior dummies have %d columns but the design has %d.', ...
                size(prior.Xd,2), size(X,2));
        end
        Yr = [Yr; prior.Yd];
        X  = [X;  prior.Xd];
    end

    % OLS, equation by equation is identical to multivariate LS here.
    Beta = (X.'*X) \ (X.'*Yr);          % (1+Kp) x K
    Uall = Yr - X*Beta;                  % residuals on the augmented sample
    U    = Uall(1:Tobs, :);              % ... and on the data rows alone

    % Degrees-of-freedom corrected covariance (note's convention).
    nreg   = size(X,2);
    % With a prior, Sigma_u is the posterior mean of the inverse-Wishart,
    % which uses the AUGMENTED residual cross-product and its own degrees of
    % freedom. Without one this reduces to the usual OLS covariance.
    if isempty(prior)
        Sigmau = (U.'*U) / (Teff - nreg);
    else
        Sigmau = (Uall.'*Uall) / (size(X,1) - nreg);
    end
    Sigmau = 0.5*(Sigmau + Sigmau.');    % enforce symmetry

    % Unpack coefficients.
    Bt = Beta.';                         % K x (1+Kp)
    if includeConst
        nu = Bt(:,1);
        Aall = Bt(:,2:end);
    else
        nu = zeros(K,1);
        Aall = Bt;
    end
    if qExo > 0
        Bexo = Aall(:, K*p+1 : K*p+qExo).';     % q x K
        Aall = Aall(:, 1:K*p);
    else
        Bexo = zeros(0, K);
    end
    Acell = cell(1,p);
    for i = 1:p
        Acell{i} = Aall(:, (i-1)*K+1 : i*K);
    end

    % Companion form  Y_t = mu + A Y_{t-1} + U_t   [Eq. (2)].
    A = zeros(K*p, K*p);
    A(1:K, :) = cell2mat(Acell);                 % top block row [A1 ... Ap]
    if p > 1
        A(K+1:end, 1:K*(p-1)) = eye(K*(p-1));    % sub-diagonal identities
    end
    J  = [eye(K), zeros(K, K*(p-1))];            % selector, y_t = J*Y_t
    mu = J.'*nu;                                 % mu = J'*nu

    model = struct('nu',nu,'Acell',{Acell},'A',A,'mu',mu,'J',J, ...
                   'Sigmau',Sigmau,'U',U,'Beta',Beta,'Bexo',Bexo,'qExo',qExo, ...
                   'K',K,'p',p,'T',T,'Teff',Teff,'includeConst',includeConst, ...
                   'prior',prior, 'XtXinv',inv(X.'*X), 'S',Uall.'*Uall, ...
                   'Tstar',size(X,1), 'nreg',size(X,2));
end
