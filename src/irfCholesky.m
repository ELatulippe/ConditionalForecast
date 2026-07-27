function IRF = irfCholesky(R, H, opts)
% IRFCHOLESKY  Structural impulse responses under the recursive identification.
%
%   IRF = irfCholesky(R, 48)
%   IRF = irfCholesky(R, 48, struct('scaleTo',1,'boot',500))
%
% Theta_j = Phi_j D with D = chol(Sigma_u,'lower'), i.e. exactly the objects
% the variance decomposition and the BK conditional forecast already use.
%
% TWO CHOICES THAT DECIDE WHETHER YOU CAN READ THE PLOT
%
% Cumulation. Variables entering in log-differences produce responses in
% GROWTH. A price puzzle is a statement about the price LEVEL, so those are
% cumulated by default: the reported path is the percent deviation of the
% level from its no-shock path. Level variables (the policy rate,
% unemployment) are left alone and read in percentage points. Set
% .cumulate = false to see the raw growth responses.
%
% Normalisation. A raw impulse response is per ONE STANDARD DEVIATION of the
% shock, and if that shock is small the response looks small for reasons that
% have nothing to do with transmission -- a smooth market rate can have a
% policy innovation of a few basis points. Set .scaleTo to rescale so the
% shocked variable moves by that amount on impact: .scaleTo = 1 with the rate
% gives responses per 100bp, which is the number to compare against the
% literature.
%
% INPUTS
%   R     output struct from main_risk_scenarios (uses .model, .data, .tr)
%   H     horizon in months
%   opts  .cumulate  cumulate dlog responses to levels (default true)
%         .scaleTo   rescale each shock so its own variable moves by this
%                    much on impact ([] = leave at one standard deviation)
%         .boot      replications for bands (0 = none, default). When the
%                    model was estimated with a Minnesota prior these are
%                    draws from the normal-inverse-Wishart POSTERIOR, not a
%                    residual bootstrap -- resampling around a shrunk point
%                    estimate would ignore the prior's own contribution to
%                    posterior uncertainty. Set .forceBoot = true to bootstrap
%                    anyway.
%         .forceBoot bootstrap even when a prior is present (default false)
%         .quantiles default [0.16 0.84] -- a 68% band
%         .seed      RNG seed (default 20260724)
%         .verbose   default true
%
% OUTPUT struct IRF
%   .resp   K x K x H   (variable, shock, horizon) in display units
%   .lo,.hi K x K x H   band, present when .boot > 0
%   .scale  K x 1       the impact move of each shock's own variable, before
%                       rescaling -- read this before judging magnitudes
%   .names, .tcode, .cumulate, .scaleTo, .horizon
%
% Bands use the same residual bootstrap as wzBands: the exogenous block is
% held fixed, residuals are resampled from months no impulse dummy switches
% on, and explosive draws are redrawn. Each draw is renormalised by its own
% impact scale, so a 100bp band really is a 100bp band.

    if nargin < 3, opts = struct(); end
    opts = withDefaults(opts, struct('cumulate',true, 'scaleTo',[], ...
        'boot',0, 'quantiles',[0.16 0.84], 'seed',20260724, 'verbose',true, ...
        'forceBoot',false));

    model = R.model;  K = model.K;  p = model.p;
    names = R.data.names;  tcode = R.data.tcode;

    [resp, scale] = oneIRF(model, H, K, tcode, opts);

    IRF = struct('resp',resp, 'scale',scale, 'names',{names}, ...
                 'tcode',{tcode}, 'cumulate',opts.cumulate, ...
                 'scaleTo',opts.scaleTo, 'horizon',H, ...
                 'quantiles',opts.quantiles);

    if opts.boot <= 0, return; end

    % ---- bootstrap bands -------------------------------------------------
    Y = R.tr.Y;  T = size(Y,1);
    Xexo = [];
    if isfield(R,'Xexo'), Xexo = R.Xexo; end
    setSeed(opts.seed);

    pool = true(model.Teff, 1);
    if ~isempty(Xexo)
        Xw = Xexo(p+1:T, :);
        isImp = sum(Xw ~= 0, 1) == 1;          % impulse dummies only
        if any(isImp), pool = ~any(Xw(:, isImp) ~= 0, 2); end
    end

    usePost = isfield(model,'prior') && ~isempty(model.prior) && ~opts.forceBoot;
    if opts.verbose
        if usePost
            fprintf('  bands from %d posterior draws (Minnesota prior).\n', opts.boot);
        else
            fprintf('  bands from %d bootstrap replications.\n', opts.boot);
        end
    end

    N = opts.boot;
    draws = zeros(K, K, H, N);
    nRej = 0;
    for d = 1:N
        for attempt = 1:200
            if usePost
                md = varPosteriorDraw(model);
            else
                Ys = simulateFrom(model, Y, p, Xexo, pool);
                md = estimateVAR(Ys, p, true, Xexo);
            end
            if max(abs(eig(md.A))) < 0.9999, break; end
            nRej = nRej + 1;
        end
        draws(:,:,:,d) = oneIRF(md, H, K, tcode, opts);
        if opts.verbose && mod(d, max(1,floor(N/5))) == 0
            fprintf('  irf bootstrap: %d/%d\n', d, N);
        end
    end

    qs = opts.quantiles;
    lo = zeros(K,K,H);  hi = zeros(K,K,H);
    for v = 1:K
        for sK = 1:K
            for h = 1:H
                q = qtile(squeeze(draws(v,sK,h,:)), qs);
                lo(v,sK,h) = q(1);  hi(v,sK,h) = q(end);
            end
        end
    end
    IRF.lo = lo;  IRF.hi = hi;  IRF.nBoot = N;  IRF.nRejected = nRej;
    IRF.bandSource = 'bootstrap';
    if usePost, IRF.bandSource = 'posterior'; end
    if opts.verbose && nRej > 0.2*(N+nRej)
        fprintf(2, ['  %.0f%% of IRF bootstrap draws were explosive and ' ...
                    'redrawn; bands are conditional on stationarity.\n'], ...
                100*nRej/(N+nRej));
    end
end

% ======================================================================
function [resp, scale] = oneIRF(model, H, K, tcode, opts)
    D  = chol(model.Sigmau, 'lower');
    MA = varMA(model, H, D);
    resp = zeros(K, K, H);
    for h = 1:H
        resp(:,:,h) = MA.Theta(:,:,h);
    end

    % Cumulate the log-difference variables into level responses (x100 = %).
    if opts.cumulate
        for v = 1:K
            if strcmp(tcode{v}, 'dlog')
                resp(v,:,:) = 100 * cumsum(resp(v,:,:), 3);
            end
        end
    end

    % Impact move of each shock's own variable, in that variable's own
    % reported units -- this is what .scaleTo divides by.
    scale = zeros(K,1);
    for sK = 1:K
        scale(sK) = resp(sK, sK, 1);
    end
    if ~isempty(opts.scaleTo)
        for sK = 1:K
            if abs(scale(sK)) > 0
                resp(:,sK,:) = resp(:,sK,:) * (opts.scaleTo / scale(sK));
            end
        end
    end
end

% ======================================================================
function Ys = simulateFrom(m, Y, p, Xexo, pool)
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

function q = qtile(x, ps)
    x = sort(x(:));  n = numel(x);
    if n == 1, q = repmat(x, 1, numel(ps)); return; end
    grid = ((1:n) - 0.5) / n;
    q = zeros(1, numel(ps));
    for j = 1:numel(ps)
        pj = ps(j);
        if pj <= grid(1),       q(j) = x(1);
        elseif pj >= grid(end), q(j) = x(end);
        else
            i = find(grid <= pj, 1, 'last');
            w = (pj - grid(i)) / (grid(i+1) - grid(i));
            q(j) = (1-w)*x(i) + w*x(i+1);
        end
    end
end

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
