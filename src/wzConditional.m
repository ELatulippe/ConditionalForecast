function [Fcond, uhat, info] = wzConditional(model, Fbase, cvars, cpaths, H, MA, mode, csd)
% WZCONDITIONAL  Waggoner-Zha (1999) conditional forecast.
%
%   [Fcond, uhat, info] = wzConditional(model, Fbase, cvars, cpaths, H, MA)
%   [...] = wzConditional(model, Fbase, cvars, cpaths, H, MA, 'reduced')
%   [...] = wzConditional(model, Fbase, cvars, cpaths, H, MA, [], csd)  % soft
%
% Impose a path for one or more variables and find the most likely shock
% sequence consistent with it, then propagate that sequence to every variable.
%
% WHICH NORM IS MINIMISED
%   The constraint R s = r has infinitely many solutions, and the scenario
%   selects one by minimising a norm. Doan, Litterman and Sims (1984) -- whom
%   the paper invokes for this step -- minimise the sum of squares of the
%   STRUCTURAL shocks, which have identity covariance. That is the "most
%   likely combination of shocks" reading: under u ~ N(0, Sigma_u) the most
%   likely u satisfying the constraint minimises the Mahalanobis norm
%   u' Sigma_u^{-1} u, equivalently eps'eps with u = D eps.
%
%   Minimising the plain Euclidean norm of the REDUCED-FORM innovations
%   instead treats every equation as equally costly to shock. With Canadian
%   monthly data the innovation standard deviations span two orders of
%   magnitude -- oil around 0.09 in logs against CPI around 0.003 -- so the
%   unweighted solution buys cheap-looking reductions in ||u|| by loading
%   enormous shocks, in standard-deviation terms, onto the low-variance
%   equations. The constrained variable still follows its path exactly, so
%   the error is invisible in the panel you constrained and shows up as
%   implausible paths for prices and housing.
%
%   Both solutions satisfy R s = r exactly. They differ in every unconstrained
%   variable. 'structural' (the default) is the correct one; 'reduced'
%   reproduces the earlier behaviour for comparison.
%
% INPUTS
%   model    struct from estimateVAR
%   Fbase    H x K baseline forecast (from baselineForecast)
%   cvars    1 x M vector of constrained variable indices (columns of y_t)
%   cpaths   H x M matrix of TARGET values (transformed units) for those
%            variables over horizons 1..H. NaN leaves a cell unconstrained.
%   H        horizon
%   MA       (optional) varMA(model,H,D) output. Supplying one built WITH D
%            avoids recomputing the structural map.
%   mode     'structural' (default) or 'reduced'
%   csd      optional, same shape as cpaths: the standard deviation of a SOFT
%           condition, in the same units as the path. 0 or NaN means hard.
%
% SOFT CONDITIONS
%   Waggoner and Zha allow conditioning on a DISTRIBUTION for the constrained
%   variable rather than a point path -- "oil around 120, give or take" rather
%   than exactly 120. Write the restriction as an observation with noise,
%
%       r = R eps + e,     e ~ N(0, Omega),   Omega = diag(csd.^2)
%
%   and update the prior eps ~ N(0, I) in the usual Gaussian way:
%
%       eps | r ~ N( R'(RR' + Omega)^{-1} r ,  I - R'(RR' + Omega)^{-1} R )
%
%   Omega = 0 recovers the hard condition exactly, with mean pinv(R) r and
%   covariance I - pinv(R) R, so the two cases share one code path. The
%   difference that matters downstream: under a soft condition the constrained
%   variable's own band is no longer degenerate, which is more honest for
%   something labelled a risk scenario.
%
%   Units. csd is in the units of the path, so for the policy rate or
%   unemployment it is percentage points -- 'rate at 6% with sd 0.5' is
%   csd = 0.5. For a log-difference variable it is in log-growth units.
%
% OUTPUTS
%   Fcond    H x K conditional forecast (transformed units)
%   uhat     KH x 1 stacked reduced-form innovations implied by the scenario
%   info     struct with R, r, rows, mode, rank, cond, the K x H structural
%            shock matrix .eps, its norm .normEps, and .constraintResidual

    K = model.K;
    if nargin < 7 || isempty(mode), mode = 'structural'; end
    if nargin < 8, csd = []; end

    % Structural map needs D; build it if the caller did not supply one.
    D = chol(model.Sigmau, 'lower');
    if nargin < 6 || isempty(MA) || ~isfield(MA,'Mstr')
        MA = varMA(model, H, D);
    end

    switch lower(mode)
        case 'structural', M = MA.Mstr;    % dy = Mstr * eps
        case 'reduced',    M = MA.Mred;    % dy = Mred * u   (legacy)
        otherwise
            error('wzConditional:mode', ...
                'mode must be ''structural'' or ''reduced'', got "%s".', mode);
    end

    % ---- Build the stacked restriction  R s = r  (Eq. 5') ----------------
    rowsel = [];   rvec = [];   svec = [];
    for h = 1:H
        for m = 1:numel(cvars)
            tgt = cpaths(h, m);
            if ~isnan(tgt)
                v   = cvars(m);
                row = (h-1)*K + v;                    % row of the stacked dy
                rowsel(end+1,1) = row;                %#ok<AGROW>
                rvec(end+1,1)   = tgt - Fbase(h, v);  %#ok<AGROW>
                sd = 0;
                if ~isempty(csd) && ~isnan(csd(h,m)), sd = csd(h,m); end
                svec(end+1,1)   = sd;                 %#ok<AGROW>
            end
        end
    end
    R = M(rowsel, :);
    r = rvec;
    isSoft = any(svec > 0);

    % ---- Minimum-norm solution, via a rank-revealing SVD -----------------
    % shat = pinv(R) r, equal to R'(R R')^{-1} r at full row rank but stable
    % when R is ill conditioned (long horizons, near-deterministic series).
    [Usvd, Ssvd, Vsvd] = svd(R, 'econ');
    sv   = diag(Ssvd);
    tol  = max(size(R)) * eps(max(sv));
    keep = sv > tol;
    rankR = sum(keep);
    condR = sv(1) / sv(find(keep,1,'last'));

    if isSoft
        % Gaussian update; Omega = 0 would reproduce the pinv solution below.
        Omega = diag(svec.^2);
        Sfull = R*R.' + Omega;
        gain  = @(w) R.' * (Sfull \ w);          % R'(RR'+Omega)^{-1} w
        shat  = gain(r);
    else
        gain  = [];
        shat  = Vsvd(:,keep) * ((Usvd(:,keep).' * r) ./ sv(keep));
    end
    if rankR < size(R,1)
        warning('wzConditional:rankDeficient', ...
            ['The constraint system has rank %d < %d rows: the scenario is ' ...
             'not fully attainable and the shortfall is being projected away. ' ...
             'Check that the constrained path is expressed in the VAR''s ' ...
             'transformed units.'], rankR, size(R,1));
    elseif condR > 1e10
        warning('wzConditional:illConditioned', ...
            ['cond(R) = %.1e. The implied shocks are only weakly pinned down; ' ...
             'the unconstrained forecasts are sensitive to small data changes.'], condR);
    end

    % ---- Map between structural and reduced-form innovations -------------
    if strcmpi(mode, 'structural')
        Eps  = reshape(shat, K, H);          % column h = eps_{t+h}
        uhat = reshape(D * Eps, K*H, 1);     % u_{t+h} = D eps_{t+h}
    else
        uhat = shat;
        Eps  = D \ reshape(shat, K, H);      % eps_{t+h} = D^{-1} u_{t+h}
    end

    % ---- Propagate to all variables --------------------------------------
    dfull = M * shat;                        % KH x 1 stacked path deviation
    Fcond = Fbase + reshape(dfull, K, H).';  % H x K

    if isSoft
        resid = NaN;                 % a soft condition is not meant to bind
    else
        resid = norm(R*shat - r);
    end
    info = struct('R',R, 'r',r, 'rows',rowsel, 'mode',lower(mode), ...
                  'rank',rankR, 'cond',condR, 'eps',Eps, ...
                  'normEps',norm(Eps(:)), 'D',D, 'soft',isSoft, ...
                  'sd',svec, 'gain',gain, ...
                  'constraintResidual', resid);
end
