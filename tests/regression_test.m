function regression_test()
%REGRESSION_TEST  Pin the numerical behaviour of the estimation pipeline.
%
%   >> regression_test
%
% Runs the model once, OFFLINE, from the shipped panel cache (cache/panel_k9.mat)
% and checks two kinds of things:
%
%   (A) ALGEBRAIC INVARIANTS that must hold for any data, so they catch a broken
%       refactor without needing platform-specific "golden" numbers:
%         * Sigma_u is symmetric and positive definite
%         * D = chol(Sigma_u,'lower') is lower-triangular and D*D' = Sigma_u
%         * the impact impulse response equals D
%         * every FEVD row sums to 1 and every share is in [0,1]
%         * total forecast-error variance is non-decreasing in the horizon
%         * Waggoner-Zha hard conditioning hits its constraint (residual ~ 0)
%         * a counterfactual with r_cf = r_actual reproduces history exactly
%
%   (B) A GOLDEN SNAPSHOT of a handful of scalars (coefficient norm, a few FEVD
%       cells, a scenario's ||eps||). On the FIRST run it records them to
%       tests/golden/golden_regression.mat; on later runs it compares within a
%       relative tolerance and FAILS on drift. Delete that file to re-baseline
%       after an intended numerical change.
%
% Everything here is deterministic: offline data, no bootstrap, no bands. The
% function errors (non-zero exit intent) if any check fails, so it is usable in
% a CI hook: matlab -batch "addpath('tests'); regression_test".

    %% --- paths (tests/ lives one level under the project root) ----------
    thisDir = fileparts(mfilename('fullpath'));
    root    = fileparts(thisDir);
    addpath(root); addpath(fullfile(root,'src'));
    assert(exist('main_risk_scenarios','file')==2, ...
        'main_risk_scenarios.m not found under ./src.');

    reltol = 1e-6;   % cross-platform (MATLAB/Octave BLAS) relative tolerance
    abstol = 1e-8;
    T = struct('name',{},'pass',{},'detail',{});   % results accumulator

    %% --- one deterministic offline run ----------------------------------
    ob = scenario_config();
    ob.offline = true;    % never touch FRED; require the shipped cache
    ob.bands   = false;
    ob.figures = false;
    ob.gate    = false;
    % Two HARD conditioning scenarios with fixed targets (no dependence on the
    % forecast origin, so the test is self-contained).
    ob.scenarios = { ...
        struct('var','oil', 'type','ramp',  'to',100, 'months',3, 'name','oil ramp (test)'), ...
        struct('var','rate','type','policy','step',0.25,'peak',4,'hold',6,'floor',2, ...
               'name','rate policy (test)') };

    try
        R = main_risk_scenarios('fred', ob);
    catch err
        fprintf(2, ['\nCould not run offline. If this says the cache is missing,\n' ...
                    'build it once online with  freeze_panel.\n\n']);
        rethrow(err);
    end

    idx = R.data.idx;  K = R.model.K;  H = ob.H;
    Sig = R.model.Sigmau;

    %% ================= (A) algebraic invariants =========================
    % Sigma_u symmetric PSD
    T = add(T, 'Sigma_u symmetric', norm(Sig - Sig','fro') <= abstol*max(1,norm(Sig,'fro')), ...
            sprintf('asym = %.2e', norm(Sig - Sig','fro')));
    ev = eig((Sig+Sig')/2);
    T = add(T, 'Sigma_u positive definite', min(real(ev)) > 0, ...
            sprintf('min eig = %.3e', min(real(ev))));

    % Cholesky factor
    D = chol(Sig,'lower');
    T = add(T, 'D lower-triangular', norm(triu(D,1),'fro') <= abstol, ...
            sprintf('upper-tri mass = %.2e', norm(triu(D,1),'fro')));
    T = add(T, 'D*D^T = Sigma_u', norm(D*D.' - Sig,'fro') <= reltol*norm(Sig,'fro'), ...
            sprintf('resid = %.2e', norm(D*D.' - Sig,'fro')));
    if isfield(R,'D')
        T = add(T, 'R.D matches chol(Sigma_u)', norm(R.D - D,'fro') <= reltol*max(1,norm(D,'fro')), ...
                sprintf('diff = %.2e', norm(R.D - D,'fro')));
    end

    % Impact IRF equals D  (Theta_0 = Phi_0 D = D), raw units, no cumulation
    I0 = irfCholesky(R, H, struct('cumulate',false,'scaleTo',[],'boot',0,'verbose',false));
    resp0 = squeeze(I0.resp(:,:,1));
    T = add(T, 'impact IRF = D', norm(resp0 - D,'fro') <= reltol*max(1,norm(D,'fro')), ...
            sprintf('diff = %.2e', norm(resp0 - D,'fro')));

    % FEVD: rows sum to 1, shares in [0,1], MSPE non-decreasing
    S = R.shares;  M = R.MSPE;   % K x K x H  and  K x H
    rowSums = sum(S, 2);                       % K x 1 x H, each should be ~1
    T = add(T, 'FEVD rows sum to 1', max(abs(rowSums(:)-1)) <= 1e-8, ...
            sprintf('max |rowsum-1| = %.2e', max(abs(rowSums(:)-1))));
    T = add(T, 'FEVD shares in [0,1]', all(S(:) >= -1e-12) && all(S(:) <= 1+1e-12), ...
            sprintf('range [%.3e, %.4f]', min(S(:)), max(S(:))));
    dM = diff(M, 1, 2);
    T = add(T, 'total FEV non-decreasing in h', all(dM(:) >= -1e-8*max(1,max(M(:)))), ...
            sprintf('min step = %.2e', min(dM(:))));

    % Baseline finite
    T = add(T, 'baseline forecast finite', all(isfinite(R.Fbase(:))), '');

    %% ---- WZ hard conditioning hits its constraint ----------------------
    if isfield(R,'wzinfo')
        for k = 1:numel(R.wzinfo)
            r = R.wzinfo{k}.constraintResidual;
            T = add(T, sprintf('scenario %d hits constraint', k), r <= 1e-6, ...
                    sprintf('residual = %.2e', r));
        end
    end

    %% ---- counterfactual identity: r_cf = r_actual reproduces history ---
    Hcf = 24;
    t0  = find(R.tr.dates == datenum(2021,12,1));
    if isempty(t0), t0 = numel(R.tr.dates) - Hcf - 1; end   % fallback origin
    ixr = idx.rate;
    rAct = R.tr.Y(t0+1:t0+Hcf, ixr);
    C = wzCounterfactual(R.model, R.tr.Y, R.tr.dates, datevec_to_ym(R.tr.dates(t0)), ...
                         ixr, rAct, Hcf, [], ixr);           % identical target
    T = add(T, 'counterfactual identity (cf=actual)', ...
            norm(C.diff(:)) <= 1e-6 && C.normEps <= 1e-6, ...
            sprintf('||diff|| = %.2e, ||eps|| = %.2e', norm(C.diff(:)), C.normEps));

    %% ================= (B) golden snapshot ==============================
    g = struct();
    g.coefNorm  = norm(R.model.A, 'fro');
    g.traceSig  = trace(Sig);
    g.detD      = prod(diag(D));
    g.fevd_cpi_oil_12   = pickShare(S, idx, 'cpi',  'oil',  12);
    g.fevd_rate_rate_12 = pickShare(S, idx, 'rate', 'rate', 12);
    g.mspe_cpi_24       = M(idx.cpi, min(24,H));
    if isfield(R,'wzinfo') && ~isempty(R.wzinfo)
        g.scen1_normEps = R.wzinfo{1}.normEps;
    end

    goldDir  = fullfile(thisDir,'golden');
    goldFile = fullfile(goldDir,'golden_regression.mat');
    if exist(goldFile,'file')
        ref = load(goldFile); ref = ref.g;
        f = fieldnames(g);
        for i = 1:numel(f)
            if ~isfield(ref, f{i}), continue; end
            a = g.(f{i});  b = ref.(f{i});
            ok = abs(a-b) <= reltol*max(1,abs(b));
            T = add(T, sprintf('golden: %s', f{i}), ok, ...
                    sprintf('now %.10g vs ref %.10g', a, b));
        end
    else
        if ~exist(goldDir,'dir'), mkdir(goldDir); end
        save(goldFile, 'g');
        fprintf(2, ['\n[golden] No baseline found -- RECORDED a new one at\n    %s\n' ...
                    'Re-run regression_test to compare against it.\n'], goldFile);
    end

    %% ================= report ===========================================
    nfail = report(T);
    if nfail > 0
        error('regression_test:fail', '%d regression check(s) FAILED.', nfail);
    end
    fprintf('\nregression_test: all %d checks passed.\n', numel(T));
end

% ======================================================================
function T = add(T, name, pass, detail)
    T(end+1) = struct('name',name, 'pass',logical(pass), 'detail',detail); %#ok<AGROW>
end

function v = pickShare(S, idx, tgt, shk, h)
    v = S(idx.(tgt), idx.(shk), h);
end

function ym = datevec_to_ym(dn)
    dv = datevec(dn);  ym = [dv(1) dv(2)];
end

function nfail = report(T)
    fprintf('\n%-42s %s\n', 'check', 'result');
    fprintf('%s\n', repmat('-',1,64));
    nfail = 0;
    for i = 1:numel(T)
        if T(i).pass
            tag = 'PASS';
        else
            tag = 'FAIL'; nfail = nfail + 1;
        end
        d = T(i).detail; if ~isempty(d), d = ['   (' d ')']; end
        fprintf('%-42s %-4s%s\n', T(i).name, tag, d);
    end
    fprintf('%s\n', repmat('-',1,64));
    fprintf('%d checks, %d failed.\n', numel(T), nfail);
end
