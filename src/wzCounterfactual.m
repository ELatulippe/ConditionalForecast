function C = wzCounterfactual(model, Y, dates, origin, cvars, cpaths, H, MA, onlyShock)
% WZCOUNTERFACTUAL  Historical counterfactual by Waggoner-Zha conditioning.
%
%   C = wzCounterfactual(model, tr.Y, tr.dates, [2021 12], idx.rate, path, 24)
%   C = wzCounterfactual(..., MA, idx.rate)     % policy shocks only
%
% "What would have happened over 2022-23 if the Bank had held the rate at
% 0.25%?" Nothing in the Waggoner-Zha algebra restricts it to the forecast
% horizon; the only change is what you perturb.
%
% THE OBJECT. Over a window the data actually realised, write the outcome as
% the baseline forecast from the origin plus the realised innovations:
%
%     y_actual = Fbase + M u_actual,   so   R u_actual = r_actual
%
% where r_actual is the realised deviation of the constrained variable from
% the baseline. A counterfactual asks for a different r. The natural answer
% keeps every realised shock that the counterfactual does not require you to
% change, i.e. solves
%
%     min ||delta||  s.t.  R (u_actual + delta) = r_cf   =>   delta = pinv(R)(r_cf - r_actual)
%
% and therefore
%
%     y_cf = y_actual + M pinv(R) (r_cf - r_actual).
%
% Note what drops out: the baseline never appears. The counterfactual is a
% perturbation of HISTORY, not of a forecast, so it inherits every realised
% shock outside the span of the constraint. That is the whole point -- the
% 2022 oil shock, the supply-chain shock and so on are all still there.
%
% As everywhere else in this package the norm is taken over the STRUCTURAL
% shocks, so the perturbation is the most likely one rather than the smallest
% in an arbitrary metric.
%
% THE CAVEAT, which belongs in any write-up using this. The reduced-form VAR
% is assumed invariant to the change in policy being contemplated. That is
% precisely what the Lucas critique denies. The exercise is informative about
% the model's internal propagation, not about what the economy would truly
% have done; the more the counterfactual departs from historical experience,
% the weaker the claim. ||delta|| relative to sqrt(K*H) is a usable gauge.
%
% INPUTS
%   model    struct from estimateVAR
%   Y        T x K transformed data (tr.Y)
%   dates    T x 1 datenums aligned with Y (tr.dates)
%   origin   [year month] of the last month BEFORE the counterfactual window,
%            or a row index into Y
%   cvars    1 x M constrained variable indices
%   cpaths   H x M counterfactual targets in TRANSFORMED units; NaN frees a cell
%   H        length of the window in months
%   MA       optional varMA(model,H,D) output
%   onlyShock  optional vector of shock indices the perturbation is allowed
%           to use. See WHICH SHOCKS below.
%
% WHICH SHOCKS DO THE WORK
%   Left unrestricted, the minimum-norm perturbation spreads across ALL
%   structural shocks -- it finds the most likely way for the constrained
%   variable to have taken a different path, not the most likely POLICY
%   reason. For a rate that is ordered late in the recursive scheme, the
%   cheapest explanation for "the rate stayed low" is often "the economy was
%   weaker", so the solution loads on demand shocks and the counterfactual
%   comes back with LOWER prices, not higher. That is a coherent answer to a
%   different question, and it is the Waggoner-Zha versus Baumeister-Kilian
%   distinction the paper raises in Section 6.
%
%   Pass onlyShock = idx.rate to confine the perturbation to the monetary
%   policy shock. The counterfactual then reads "the same history, but the
%   policy shocks were whatever they had to be to hold the rate flat", which
%   is what "if the Bank had not tightened" usually means. With H constraints
%   and H free policy shocks the system is square, so the path is still hit
%   exactly; check C.normEps, which will be larger, since the perturbation
%   can no longer use the cheapest shocks available.
%
% OUTPUT struct C
%   .dates     H x 1 datenums of the window
%   .actual    H x K realised data over the window (transformed units)
%   .cf        H x K counterfactual
%   .diff      H x K  cf - actual
%   .Fbase     H x K baseline forecast from the origin, for reference
%   .eps       K x H structural shock perturbation
%   .epsNorm   K x 1 norm of the perturbation per shock -- read this to see
%              which shocks the counterfactual actually used
%   .normEps   ||delta||, and .ratio = ||delta|| / sqrt(K*H)
%   .rank, .cond, .constraintResidual

    K = model.K;  p = model.p;
    if nargin < 8 || isempty(MA) || ~isfield(MA,'Mstr')
        D  = chol(model.Sigmau, 'lower');
        MA = varMA(model, H, D);
    end
    M = MA.Mstr;

    % ---- locate the origin ---------------------------------------------
    if numel(origin) == 2
        t0 = find(dates == datenum(origin(1), origin(2), 1), 1);
        if isempty(t0)
            error('wzCounterfactual:origin', ...
                'Origin %04d-%02d is not in the sample (%s .. %s).', ...
                origin(1), origin(2), datestr(dates(1),'yyyy-mm'), ...
                datestr(dates(end),'yyyy-mm'));
        end
    else
        t0 = origin;
    end
    if t0 < p
        error('wzCounterfactual:origin', ...
            'Origin must leave %d lags before it; it is at row %d.', p, t0);
    end
    if t0 + H > size(Y,1)
        error('wzCounterfactual:window', ...
            ['The window runs %d months past the end of the sample. Shorten H ' ...
             'or move the origin back.'], t0 + H - size(Y,1));
    end

    actual = Y(t0+1 : t0+H, :);
    Fbase  = baselineForecast(model, Y(t0-p+1 : t0, :), H);

    % ---- restriction rows, and the realised value on them ---------------
    rowsel = [];  rcf = [];  ract = [];
    for h = 1:H
        for m = 1:numel(cvars)
            tgt = cpaths(h, m);
            if ~isnan(tgt)
                v = cvars(m);
                rowsel(end+1,1) = (h-1)*K + v;                 %#ok<AGROW>
                rcf(end+1,1)    = tgt          - Fbase(h, v);  %#ok<AGROW>
                ract(end+1,1)   = actual(h, v) - Fbase(h, v);  %#ok<AGROW>
            end
        end
    end
    R = M(rowsel, :);

    % ---- optionally confine the perturbation to selected shocks ---------
    if nargin >= 9 && ~isempty(onlyShock)
        colsel = [];
        for h = 1:H
            for q = 1:numel(onlyShock)
                colsel(end+1,1) = (h-1)*K + onlyShock(q);   %#ok<AGROW>
            end
        end
    else
        colsel = (1:size(R,2)).';
    end
    Rsub = R(:, colsel);

    % ---- minimal structural perturbation --------------------------------
    [Usvd, Ssvd, Vsvd] = svd(Rsub, 'econ');
    sv   = diag(Ssvd);
    tol  = max(size(Rsub)) * eps(max(sv));
    keep = sv > tol;
    dsub = Vsvd(:,keep) * ((Usvd(:,keep).' * (rcf - ract)) ./ sv(keep));
    delta = zeros(size(R,2), 1);
    delta(colsel) = dsub;

    cf = actual + reshape(M*delta, K, H).';

    C = struct();
    C.dates   = dates(t0+1 : t0+H);
    C.actual  = actual;
    C.cf      = cf;
    C.diff    = cf - actual;
    C.Fbase   = Fbase;
    C.eps     = reshape(delta, K, H);
    C.normEps = norm(delta);
    C.ratio   = C.normEps / sqrt(K*H);
    C.rank    = sum(keep);
    C.cond    = sv(1) / sv(find(keep,1,'last'));
    C.constraintResidual = norm(R*delta - (rcf - ract));
    C.rows    = rowsel;
    C.epsNorm = sqrt(sum(C.eps.^2, 2));
    C.onlyShock = [];
    if nargin >= 9, C.onlyShock = onlyShock; end
end
