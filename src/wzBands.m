function B = wzBands(Y, p, Xexo, scen, H, opts)
% WZBANDS  Uncertainty bands for Waggoner-Zha conditional forecasts.
%
%   B = wzBands(Y, p, Xexo, scen, H)
%   B = wzBands(Y, p, Xexo, scen, H, struct('nDraws',2000,'parameter',false))
%
% The paper reports conditional means only. Waggoner and Zha (1999) is largely
% about the distribution around them, so this implements the part the paper
% cites but does not use. Two sources of uncertainty are combined.
%
% SHOCK UNCERTAINTY (parameters held fixed).
%   Stack the structural shocks eps ~ N(0, I_KH), with the path deviation
%   dy = M eps for M = Mstr and the scenario written as R eps = r.
%   Conditioning a Gaussian on an affine restriction gives another Gaussian
%   supported on that subspace:
%
%       eps | R eps = r  ~  N( ehat, P ),   ehat = pinv(R) r,  P = I - pinv(R) R
%
%   ehat is exactly the point solution wzConditional returns, and P is the
%   orthogonal projector onto null(R). A draw is ehat + P z with z ~ N(0,I).
%   P is never formed: it is applied as z - pinv(R)(R z) through the same SVD
%   factors used for the point solution. Since R P = R - R pinv(R) R = 0, every
%   draw satisfies the constraint EXACTLY, so the constrained variable's band
%   has zero width at every horizon. test_engines checks that.
%
% PARAMETER UNCERTAINTY (residual bootstrap).
%   Each draw resamples the reduced-form residuals, regenerates the sample,
%   and re-estimates. R and r are therefore draw-specific -- A, Sigma_u, D,
%   the MA coefficients and the baseline all move -- so the constraint is
%   re-imposed inside the loop. Computing one ehat and adding parameter noise
%   afterwards would be wrong.
%
%   With a near-unit root (max|eig(A)| about 0.98 here) the OLS bootstrap is
%   visibly biased, so the DGP is bias-adjusted following Kilian (1998):
%   the bias is estimated by a first bootstrap pass and subtracted, and if the
%   adjusted companion is explosive the adjustment is shrunk until it is not.
%   Explosive draws in the main pass are rejected and redrawn.
%
%   When the model carries pandemic dummies, residuals are resampled from the
%   NON-dummied months only. Those months' residuals are zero by construction,
%   so including them would reintroduce the 2020 observations through the back
%   door and understate the innovation variance.
%
% INPUTS
%   Y      T x K transformed data (the VAR's input)
%   p      lag order
%   Xexo   T x q exogenous block (may be empty)
%   scen   struct array from buildScenarios (.name, .cvar, .path)
%   H      horizon
%   opts   struct, all optional:
%     .nDraws     number of draws (default 1000)
%     .quantiles  probabilities to report (default [0.05 0.16 0.50 0.84 0.95])
%     .shock      include shock uncertainty (default true)
%     .parameter  include parameter uncertainty (default true)
%     .prior      Minnesota prior struct (.Yd/.Xd) to estimate with. When one
%                 is supplied, parameter uncertainty is drawn from the
%                 normal-inverse-Wishart POSTERIOR rather than bootstrapped:
%                 resampling residuals around a shrunk point estimate would
%                 ignore the prior's own contribution to posterior
%                 uncertainty. Kilian's bias adjustment is then skipped -- it
%                 corrects OLS small-sample bias, which shrinkage has already
%                 displaced.
%     .forceBoot  bootstrap even when a prior is present (default false); the
%                 bootstrap then re-estimates WITH the prior each draw
%     .biasCorrect  Kilian adjustment of the bootstrap DGP (default true)
%     .nBias      draws for the bias estimate (default 200)
%     .seed       RNG seed (default 20260723)
%     .verbose    progress reporting (default true)
%     .keepDraws  return the raw draws as well as the quantiles (default
%                 true). Needed for plotting: the display transforms are
%                 nonlinear -- dispSeries reconstructs levels as
%                 exp(cumsum(.)) and year-over-year growth from those levels
%                 -- so each draw must be transformed and the quantiles taken
%                 afterwards, never the reverse. Costs H*K*nS*N doubles,
%                 about 11 MB at the defaults.
%
% OUTPUT struct B
%   .quantiles              the probabilities reported
%   .base       H x K x nQ  bands for the unconditional forecast
%   .scen(k).name
%   .scen(k).cond  H x K x nQ   bands for the conditional forecast
%   .scen(k).diff  H x K x nQ   bands for (conditional mean - baseline), which
%                               carry PARAMETER uncertainty only: with the
%                               parameters fixed that difference is
%                               deterministic. This is the object to read when
%                               asking whether a scenario is distinguishable
%                               from the baseline, because the shock
%                               uncertainty common to both cancels.
%   .condDraws  H x K x nS x N   conditional draws, shock noise included
%   .baseDraws  H x K x N        baseline draws, shock noise included
%   .condMean   H x K x nS x N   per-draw conditional MEAN (no shock noise)
%   .baseMean   H x K x N        per-draw baseline MEAN (no shock noise)
%
%   The last two exist because of a trap in reading the fan chart. The point
%   baseline computed from the OLS estimate is NOT the centre of the
%   bootstrap distribution: the DGP is bias-adjusted and explosive draws are
%   rejected, both of which shift the baseline across draws. Comparing the OLS
%   baseline with the bootstrap median of the conditional forecast therefore
%   mixes the scenario effect with that shift. Compare .condMean against
%   .baseMean draw by draw instead -- the shift is common to both and cancels.
%   .nDraws, .nRejected, .biasShrink

    [T, K] = size(Y);
    if nargin < 6, opts = struct(); end
    opts = withDefaults(opts, struct('nDraws',1000, ...
        'quantiles',[0.05 0.16 0.50 0.84 0.95], 'shock',true, ...
        'parameter',true, 'biasCorrect',true, 'nBias',200, ...
        'seed',20260723, 'verbose',true, 'keepDraws',true, ...
        'prior',[], 'forceBoot',false));
    setSeed(opts.seed);

    qs = opts.quantiles(:).';  nQ = numel(qs);  nS = numel(scen);  N = opts.nDraws;

    model0  = estimateVAR(Y, p, true, Xexo, opts.prior);
    usePost = ~isempty(opts.prior) && ~opts.forceBoot;
    Ylast  = Y(end-p+1:end, :);

    % Rows usable for residual resampling: exclude any month the exogenous
    % block switches on, whose residual is zero by construction.
    % Rows usable for resampling: exclude months switched on by an IMPULSE
    % dummy (a column with a single nonzero, e.g. a pandemic month), whose
    % residual is zero by construction. Recurring regressors such as month
    % dummies are active in most rows and must NOT empty the pool.
    pool = true(model0.Teff, 1);
    if ~isempty(Xexo)
        Xw    = Xexo(p+1:T, :);
        isImp = sum(Xw ~= 0, 1) == 1;
        if any(isImp)
            pool = ~any(Xw(:, isImp) ~= 0, 2);
        end
    end
    if sum(pool) < 20
        error('wzBands:pool', 'Only %d residuals available to resample.', sum(pool));
    end

    % ---- bias-adjusted DGP for the bootstrap ---------------------------
    dgp = model0;  shrink = NaN;
    if opts.parameter && opts.biasCorrect && ~usePost
        [dgp, shrink] = biasAdjust(model0, Y, p, Xexo, pool, opts);
        if opts.verbose
            fprintf('  Kilian bias adjustment applied (shrink factor %.2f).\n', shrink);
        end
    end

    condDraw = zeros(H, K, nS, N);
    diffDraw = zeros(H, K, nS, N);
    baseDraw = zeros(H, K, N);
    condMean = zeros(H, K, nS, N);   % Fb_d + M_d ehat_d, no shock noise
    baseMean = zeros(H, K, N);       % Fb_d itself
    nRej = 0;

    for d = 1:N
        if opts.parameter
            if usePost
                [md, rej] = drawPosterior(model0);
            else
                [md, rej] = drawModel(dgp, Y, p, Xexo, pool, opts.prior);
            end
            nRej = nRej + rej;
        else
            md = model0;
        end

        Dd  = chol(md.Sigmau, 'lower');
        MAd = varMA(md, H, Dd);
        M   = MAd.Mstr;
        Fb  = baselineForecast(md, Ylast, H);

        baseMean(:,:,d) = Fb;
        zb = zeros(K*H,1);
        if opts.shock, zb = randn(K*H,1); end
        baseDraw(:,:,d) = Fb + reshape(M*zb, K, H).';

        for k = 1:nS
            [R, r, svec] = buildRestriction(M, Fb, scen(k), H, K);
            soft = any(svec > 0);
            if soft
                % Gaussian update. A draw is ehat + (z - G(Rz + e)) with
                % z ~ N(0,I) and e ~ N(0,Omega); its covariance is exactly
                % I - R'(RR'+Omega)^{-1} R, so no factorisation is needed and
                % Omega = 0 collapses to the projector used below.
                Sfull = R*R.' + diag(svec.^2);
                G     = @(w) R.' * (Sfull \ w);
                ehat  = G(r);
            else
                [ehat, applyPinv] = minNorm(R, r);
            end
            diffDraw(:,:,k,d) = reshape(M*ehat, K, H).';
            condMean(:,:,k,d) = Fb + diffDraw(:,:,k,d);
            if opts.shock
                z = randn(K*H,1);
                if soft
                    e = ehat + (z - G(R*z + svec.*randn(numel(svec),1)));
                else
                    e = ehat + (z - applyPinv(R*z));   % project onto null(R)
                end
            else
                e = ehat;
            end
            condDraw(:,:,k,d) = Fb + reshape(M*e, K, H).';
        end

        if opts.verbose && mod(d, max(1,floor(N/10))) == 0
            fprintf('  bands: %d/%d draws\n', d, N);
        end
    end

    B = struct('quantiles',qs, 'nDraws',N, 'nRejected',nRej, 'biasShrink',shrink);
    B.bandSource = 'bootstrap';
    if usePost, B.bandSource = 'posterior'; end
    rejRate = nRej / max(N + nRej, 1);
    if opts.parameter && ~usePost && rejRate > 0.20
        warning('wzBands:manyExplosive', ...
            ['%.0f%% of bootstrap replications came out explosive and were ' ...
             'redrawn. With a near-unit root the bias adjustment pushes the ' ...
             'DGP toward the unit circle, so a large share of draws is ' ...
             'discarded and the bands are conditional on stationarity -- ' ...
             'read them as understating parameter uncertainty at long ' ...
             'horizons. Setting biasCorrect = false lowers the rejection ' ...
             'rate at the cost of a downward-biased persistence estimate.'], ...
             100*rejRate);
    end
    B.base = quantAlong(baseDraw, qs);
    B.scen = struct('name',{},'cond',{},'diff',{});
    for k = 1:nS
        B.scen(k).name = scen(k).name;
        B.scen(k).cond = quantAlong(squeeze(condDraw(:,:,k,:)), qs);
        B.scen(k).diff = quantAlong(squeeze(diffDraw(:,:,k,:)), qs);
    end
    if opts.keepDraws
        B.condDraws = condDraw;
        B.baseDraws = baseDraw;
        B.condMean  = condMean;
        B.baseMean  = baseMean;
    end
end

% ======================================================================
function [R, r, svec] = buildRestriction(M, Fbase, sc, H, K)
    cvars = sc.cvar;  cpaths = sc.path;
    hasSd = isfield(sc,'sd') && ~isempty(sc.sd);
    rowsel = [];  rvec = [];  svec = [];
    for h = 1:H
        for m = 1:numel(cvars)
            tgt = cpaths(h, m);
            if ~isnan(tgt)
                rowsel(end+1,1) = (h-1)*K + cvars(m);          %#ok<AGROW>
                rvec(end+1,1)   = tgt - Fbase(h, cvars(m));    %#ok<AGROW>
                sd = 0;
                if hasSd && ~isnan(sc.sd(h,m)), sd = sc.sd(h,m); end
                svec(end+1,1)   = sd;                          %#ok<AGROW>
            end
        end
    end
    R = M(rowsel, :);  r = rvec;
end

function [ehat, applyPinv] = minNorm(R, r)
% Minimum-norm solution and a closure applying pinv(R) to any vector, sharing
% the SVD so the projector never has to be formed.
    [U, S, V] = svd(R, 'econ');
    sv = diag(S);
    tol = max(size(R)) * eps(max(sv));
    keep = sv > tol;
    Uk = U(:,keep);  Vk = V(:,keep);  svk = sv(keep);
    ehat = Vk * ((Uk.' * r) ./ svk);
    applyPinv = @(w) Vk * ((Uk.' * w) ./ svk);
end

% ======================================================================
function [md, rej] = drawPosterior(model0)
% One draw from the normal-inverse-Wishart posterior, rejecting explosive ones.
    rej = 0;
    for attempt = 1:200
        md = varPosteriorDraw(model0);
        if max(abs(eig(md.A))) < 0.9999, return; end
        rej = rej + 1;
    end
    error('wzBands:explosive', ...
        'Could not draw a stationary posterior sample in 200 attempts.');
end

function [md, rej] = drawModel(dgp, Y, p, Xexo, pool, prior)
% One residual-bootstrap replication, rejecting explosive draws. The prior, if
% any, is reapplied to each draw so the bootstrap resamples the same estimator.
    if nargin < 6, prior = []; end
    rej = 0;
    for attempt = 1:200
        Ys = simulateFrom(dgp, Y, p, Xexo, pool);
        md = estimateVAR(Ys, p, true, Xexo, prior);
        if max(abs(eig(md.A))) < 0.9999, return; end
        rej = rej + 1;
    end
    error('wzBands:explosive', ...
        'Could not draw a stationary bootstrap replication in 200 attempts.');
end

function Ys = simulateFrom(m, Y, p, Xexo, pool)
% Regenerate a sample from m, keeping the first p observations and the
% exogenous block fixed, resampling residuals from the allowed rows only.
    [T, K] = size(Y);
    idxPool = find(pool);
    pick = idxPool(ceil(rand(T-p,1)*numel(idxPool)));
    Us = m.U(pick, :);

    Ys = zeros(T, K);
    Ys(1:p, :) = Y(1:p, :);
    for t = p+1:T
        yt = m.nu;
        for i = 1:p
            yt = yt + m.Acell{i} * Ys(t-i, :).';
        end
        if m.qExo > 0
            yt = yt + m.Bexo.' * Xexo(t, :).';
        end
        Ys(t, :) = (yt + Us(t-p, :).').';
    end
end

% ======================================================================
function [dgp, delta] = biasAdjust(model0, Y, p, Xexo, pool, opts)
% Kilian (1998): estimate the bootstrap bias in the autoregressive
% coefficients, subtract it, and shrink the adjustment if that makes the
% companion explosive.
    K = model0.K;
    acc = zeros(K, K*p);  nOK = 0;
    for b = 1:opts.nBias
        Ys = simulateFrom(model0, Y, p, Xexo, pool);
        mb = estimateVAR(Ys, p, true, Xexo);
        acc = acc + cell2mat(mb.Acell);  nOK = nOK + 1;
    end
    Abar = acc / nOK;
    Aols = cell2mat(model0.Acell);
    bias = Abar - Aols;

    delta = 1;
    while delta > 0
        Aadj = Aols - delta*bias;
        if max(abs(eig(companionOf(Aadj, K, p)))) < 0.9999, break; end
        delta = delta - 0.01;
    end
    if delta <= 0, delta = 0; Aadj = Aols; end

    dgp = model0;
    for i = 1:p
        dgp.Acell{i} = Aadj(:, (i-1)*K+1 : i*K);
    end
    dgp.A = companionOf(Aadj, K, p);
end

function A = companionOf(Atop, K, p)
    A = zeros(K*p, K*p);
    A(1:K, :) = Atop;
    if p > 1, A(K+1:end, 1:K*(p-1)) = eye(K*(p-1)); end
end

% ======================================================================
function Q = quantAlong(X, qs)
% X is H x K x N; return H x K x numel(qs).
    [H, K, N] = size(X);
    Q = zeros(H, K, numel(qs));
    for h = 1:H
        for k = 1:K
            Q(h,k,:) = qtile(squeeze(X(h,k,:)), qs);
        end
    end
end

function q = qtile(x, ps)
% Linear-interpolation quantiles on the (i-0.5)/n grid, matching the default
% MATLAB definition. Written out so no toolbox is needed.
    x = sort(x(:));  n = numel(x);
    if n == 1, q = repmat(x, 1, numel(ps)); return; end
    grid = ((1:n) - 0.5) / n;
    q = zeros(1, numel(ps));
    for j = 1:numel(ps)
        pj = ps(j);
        if pj <= grid(1)
            q(j) = x(1);
        elseif pj >= grid(end)
            q(j) = x(end);
        else
            i = find(grid <= pj, 1, 'last');
            w = (pj - grid(i)) / (grid(i+1) - grid(i));
            q(j) = (1-w)*x(i) + w*x(i+1);
        end
    end
end

% ======================================================================
function setSeed(s)
    if exist('rng','builtin') || exist('rng','file')
        try rng(s); return; catch, end
    end
    randn('seed', s);  rand('seed', s);
end

function o = withDefaults(o, d)
    f = fieldnames(d);
    for i = 1:numel(f)
        if ~isfield(o, f{i}), o.(f{i}) = d.(f{i}); end
    end
end
