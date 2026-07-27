function R = main_risk_scenarios(source, opts)
% MAIN_RISK_SCENARIOS  Replicate Moran-Stevanovic-Surprenant (2025) risk
% scenarios for the Canadian economy, following the derivation note.
%
%   R = main_risk_scenarios('fred', struct('houseFile','nhpi_monthly.csv'))
%   R = main_risk_scenarios('local', struct('csvDir','./csv', ...
%                                           'houseFile','nhpi_monthly.csv'))
%   R = main_risk_scenarios('synthetic');   % offline demo with fake data
%
% The run proceeds in the order in which things can go wrong:
%
%   1. load the panel and ASSERT the terminal month. If the sample was
%      silently trimmed, every scenario would start from the wrong row, so
%      this is checked before anything is estimated.
%   2. estimate the VAR at the paper's lag order (p = 3, fixed; AIC is
%      reported alongside as a diagnostic only).
%   3. compute and PRINT the forecast-error variance decomposition, and
%      compare it with Table 1. Table 1 depends only on the estimated VAR,
%      so it is the cleanest test of the data layer. If it fails, the
%      scenario figures are skipped by default -- there is nothing to learn
%      from them until the VAR itself matches.
%   4. only then: the four WZ conditional forecasts, the BK robustness
%      check, and Figures 1-4.
%
% OPTS (struct, all fields optional)
%   .houseFile  monthly house-price CSV -- see loadCanadaData. Required for
%               'fred'/'local' unless .allowHouseSpline is set.
%   .cpiFile    seasonally adjusted CPI CSV (the FRED default is NOT SA);
%               StatCan table 18-10-0006-01 downloads in a supported layout
%   .cpiSeries  which row to take from a wide CPI export (default 'All-items')
%   .rateFile   monthly policy/short-rate CSV, overriding the FRED Bank Rate
%   .rateSeries which row to take from a wide rate export
%   .files      struct keyed by variable name giving a CSV for ANY series,
%               e.g. struct('unemp','1410028701-eng.csv'). Use this when a
%               FRED series stops updating and caps the panel.
%   .rows       struct keyed the same way, giving the wide-export row label
%   .extra      cell array of structs ADDING variables to the panel; see
%               loadCanadaData for the fields and for the ordering caveat
%   .ids        struct keyed by variable name overriding the FRED series ID,
%               e.g. struct('unemp','LRUN64TTCAM156S')
%   .csvDir, .cacheFile, .allowHouseSpline   passed through to loadCanadaData
%   .p          lag order (default 3, the paper's choice)
%   .H          forecast horizon in months (default 48)
%   .sample     [startY startM endY endM] (default [1992 1 2022 12], the
%               paper's window). Give just [startY startM] to run to the
%               LATEST month the data support -- a real-time conditional
%               forecast rather than a replication. The terminal-date
%               assertion is then skipped (there is nothing to assert
%               against) and the Table 1 gate is switched off, since the
%               paper's decomposition was estimated on its own sample and is
%               not a valid target for a different one.
%   .requireEnd logical, error if the panel does not reach the requested
%               terminal month (default true)
%   .prior      [] (default, OLS) or a struct passed to minnesotaDummies,
%               e.g. struct('lambda',0.2). A Minnesota prior lets you run a
%               long lag order without spending a free parameter on every
%               coefficient; with K = 9 and p = 12 that is 123 regressors per
%               equation on about 400 months, which OLS cannot support.
%               Uncertainty is then sampled from the posterior rather than
%               bootstrapped -- see irfCholesky.
%   .seasonalDummies  logical, default false. Append eleven month dummies to
%               the exogenous block, removing deterministic monthly
%               seasonality from every equation without touching the data.
%               Use it when seasonalityCheck flags a series you want to keep
%               -- the Bank of Canada non-energy commodity index, for one,
%               whose agricultural and fish components are genuinely seasonal.
%   .covid      pandemic dummies: [] (default, none), 'default' for
%               March-August 2020, or [y1 m1; y2 m2]. DIAGNOSTIC ONLY -- the
%               paper uses no such dummies. See covidDummies.
%   .gate       logical, skip the figures when the Table 1 check fails
%               (default true; set false to draw them anyway)
%   .figures    logical, draw figures at all (default true)
%   .paperFigs  logical, draw the paper-replication figures 1-4 (default
%               true). Set false to produce only the band fan charts.
%   .diffFigs   logical, draw the scenario-minus-baseline figures 6
%               (default true; requires .bands)
%   .figStart   [year month] where the figure x-axis begins (default
%               [2019 1], matching the paper's Figures 1-3)
%   .figEnd     [year month] where it ends (default: the full horizon). Use
%               this to compute at H = 48 and plot two years -- the estimate
%               is unchanged, only the view is cropped. Shortening H instead
%               is a different exercise: R loses rows, so the scenario asks
%               less of the model and ||eps|| falls.
%   .figSub     subfolder name under figures/ for this run's output
%               (default '' = figures/ itself). E.g. 'temporary_shocks'
%               writes to figures/temporary_shocks/. Created if absent.
%   .scenarios  scenario spec passed to buildScenarios (default [] = the
%               paper's four). See buildScenarios for the field reference.
%   .bands      logical, compute uncertainty bands for the conditional
%               forecasts and draw a fan chart per scenario (default false --
%               it costs a bootstrap). See wzBands.
%   .bandOpts   struct passed straight to wzBands (.nDraws, .parameter,
%               .shock, .biasCorrect, ...)
%
% Returns struct R with the model, forecasts, decomposition and Table 1 check.

    if nargin < 1 || isempty(source), source = 'fred'; end
    if nargin < 2, opts = struct(); end
    opts = withDefaults(opts, struct('p',3, 'H',48, 'sample',[1992 1 2022 12], ...
                                     'requireEnd',true, 'gate',true, 'figures',true, ...
                                     'paperFigs',true, 'diffFigs',true, ...
                                     'covid',[], 'seasonalDummies',false, ...
                                     'prior',[], ...
                                     'figStart',[2019 1], 'figEnd',[], 'figSub','', ...
                                     'bands',false, 'bandOpts',struct(), ...
                                     'scenarios',[], ...
                                     'houseFile','', 'houseSeries','', ...
                                     'cpiFile','', 'cpiSeries','', ...
                                     'rateFile','', 'rateSeries','', ...
                                     'files',struct(), 'rows',struct(), ...
                                     'ids',struct(), 'extra',{{}}, 'csvDir','', ...
                                     'cacheFile','', 'allowHouseSpline',false, ...
                                     'dataDir','', 'cacheDir','', 'offline',false));

    if exist('OCTAVE_VERSION','builtin')
        try graphics_toolkit('gnuplot'); catch, end
    end
    here = fileparts(mfilename('fullpath')); addpath(here);

    % Locate the project root (this file lives in <root>/src). Data, caches and
    % figures are resolved against it, so the scripts run from any working dir.
    if exist('projectRoot','file') == 2, root = projectRoot(); else, root = here; end

    % Default the data/ and cache/ folders when they exist and were not given.
    if isempty(opts.dataDir) && exist(fullfile(root,'data'),'dir')
        opts.dataDir = fullfile(root,'data');
    end
    if isempty(opts.cacheDir) && exist(fullfile(root,'cache'),'dir')
        opts.cacheDir = fullfile(root,'cache');
    end

    outdir = fullfile(root,'figures'); if ~exist(outdir,'dir'), mkdir(outdir); end
    if ~isempty(opts.figSub)
        outdir = fullfile(outdir, opts.figSub);
        if ~exist(outdir,'dir'), mkdir(outdir); end
    end

    H = opts.H;  smp = opts.sample;

    % [startY startM] alone means "run to the latest available month".
    openEnded = (numel(smp) == 2) || any(~isfinite(smp(3:4)));
    if openEnded
        cv  = clock;
        smp = [smp(1) smp(2) cv(1) cv(2)];      % ask for everything; the
        opts.requireEnd = false;                % loader trims to what exists
        fprintf(['Open-ended sample: requesting data through %04d-%02d and ' ...
                 'using whatever\nthe series support.\n'], smp(3), smp(4));
    end

    % ================= 1. data, and the terminal-month check =============
    fprintf('Loading data (%s)...\n', source);
    dopts = struct('houseFile',opts.houseFile, 'houseSeries',opts.houseSeries, ...
                   'cpiFile',opts.cpiFile, 'cpiSeries',opts.cpiSeries, ...
                   'rateFile',opts.rateFile, 'rateSeries',opts.rateSeries, ...
                   'files',opts.files, 'rows',opts.rows, 'ids',opts.ids, ...
                   'extra',{opts.extra}, ...
                   'openEnded',openEnded, ...
                   'csvDir',opts.csvDir, 'cacheFile',opts.cacheFile, ...
                   'dataDir',opts.dataDir, 'cacheDir',opts.cacheDir, ...
                   'offline',opts.offline, ...
                   'allowHouseSpline',opts.allowHouseSpline);
    [data, span] = loadCanadaData(smp(1),smp(2),smp(3),smp(4), source, dopts);

    % Where the cache actually lives (bare name -> cache/), for the checks below.
    cachePath = opts.cacheFile;
    if ~isempty(cachePath) && ~isempty(opts.cacheDir) && isempty(fileparts(cachePath))
        cachePath = fullfile(opts.cacheDir, cachePath);
    end

    wantEnd = datenum(smp(3), smp(4), 1);
    gotEnd  = data.dates(end);
    if ~openEnded && ~isempty(cachePath) && exist(cachePath,'file') ...
            && gotEnd < wantEnd - 20
        fprintf(2, ['\n*** The cache %s ends %s, short of the requested %s.\n' ...
                    '    cacheFile is read before the date range is applied, so this is\n' ...
                    '    the OLD panel. Delete it or use a new cacheFile name. ***\n\n'], ...
                opts.cacheFile, datestr(gotEnd,'yyyy-mm'), datestr(wantEnd,'yyyy-mm'));
    end
    fprintf('\nPanel: %s .. %s  (%d months, %d series)\n', ...
            datestr(data.dates(1),'yyyy-mm'), datestr(gotEnd,'yyyy-mm'), ...
            numel(data.dates), size(data.levels,2));
    if ~strcmpi(source,'synthetic') && ~openEnded && gotEnd ~= wantEnd
        short = {span([span.last] < wantEnd).name};
        msg = sprintf(['Panel ends %s but the paper''s sample ends %s.\n' ...
               'Short series: %s\n' ...
               'Every scenario reads its starting level off the LAST row ' ...
               '(oil0, u0, r0), so the whole exercise would be shifted.\n' ...
               'Several of the OECD-sourced FRED ids used here have been ' ...
               'discontinued; substitute a live series or shorten .sample.'], ...
               datestr(gotEnd,'yyyy-mm'), datestr(wantEnd,'yyyy-mm'), ...
               strjoin(short, ', '));
        if opts.requireEnd
            error('main_risk_scenarios:terminalDate', '%s', msg);
        else
            fprintf(2, '\n*** %s ***\n\n', msg);
        end
    end
    fprintf('Last observed levels: oil %.1f, unemp %.1f%%, rate %.2f%%\n', ...
            data.levels(end,data.idx.oil), data.levels(end,data.idx.unemp), ...
            data.levels(end,data.idx.rate));

    tr  = transformData(data);
    Y   = tr.Y;  idx = tr.idx;

    seas = seasonalityCheck(tr);

    % ================= 2. estimation at the paper's lag order ============
    p = opts.p;
    [Xexo, exoLab] = covidDummies(tr.dates, opts.covid);
    covLab = exoLab;                       % before the seasonal block is added
    if opts.seasonalDummies
        dv = datevec(tr.dates);  mo = dv(:,2);
        Sd = zeros(numel(mo), 11);
        for mm = 2:12, Sd(:, mm-1) = double(mo == mm); end
        Xexo   = [Xexo, Sd];
        exoLab = [exoLab, arrayfun(@(mm) sprintf('m%02d',mm), 2:12, ...
                                   'UniformOutput', false)];
        fprintf('Seasonal dummies active (11 month indicators).\n');
    end
    if ~isempty(covLab)
        fprintf(2, ['\nPandemic dummies active for %s..%s (%d months). This is a\n' ...
                    'DIAGNOSTIC: the paper uses none, so these results are not a\n' ...
                    'literal replication.\n'], covLab{1}, covLab{end}, numel(covLab));
    end
    [pAIC, aic] = selectLag(Y, max(6, p), Xexo);
    fprintf('\nUsing p = %d (the paper''s choice). AIC would pick p = %d.\n', p, pAIC);
    if pAIC ~= p
        fprintf(['   AIC disagreeing is expected if CPI is not seasonally ' ...
                 'adjusted -- the criterion\n   then buys lags to chase the ' ...
                 'seasonal pattern. p is fixed, not selected.\n']);
    end
    priorStruct = [];
    if ~isempty(opts.prior)
        po = opts.prior;
        if ~isstruct(po), po = struct('lambda', po); end
        [Yd, Xd, pinfo] = minnesotaDummies(Y, p, tr.tcode, size(Xexo,2), po);
        priorStruct = struct('Yd',Yd, 'Xd',Xd, 'info',pinfo);
        fprintf(['Minnesota prior: lambda = %.3g, %d dummy observations ' ...
                 '(%d data rows).\n'], pinfo.lambda, pinfo.nDummy, size(Y,1)-p);
    end
    model = estimateVAR(Y, p, true, Xexo, priorStruct);
    fprintf('VAR(%d) on %d obs; max|eig(A)| = %.3f.\n', ...
            p, model.Teff, max(abs(eig(model.A))));

    D   = chol(model.Sigmau, 'lower');
    MAd = varMA(model, H, D);          % carries both Mred and Mstr

    % Table 1 is a target only for the paper's own sample. On any other
    % window the comparison is informative but not a pass/fail criterion.
    paperWindow = (gotEnd == datenum(2022,12,1)) && (data.dates(1) == datenum(1992,1,1));
    if opts.gate && ~paperWindow
        opts.gate = false;
        fprintf(['\nSample is not the paper''s 1992-01..2022-12 window, so the ' ...
                 'Table 1 gate is off.\nThe comparison below is still printed ' ...
                 'as a diagnostic.\n']);
    end

    % ================= 3. variance decomposition, then the gate ==========
    hzn = [3 12 24 48];  hzn = hzn(hzn <= H);
    [shares, MSPE] = varianceDecomp(model, H, D, hzn, ...
                                    data.names, data.names);
    chk = checkTable1(shares, idx);

    R = struct('data',data,'span',span,'tr',tr,'model',model,'p',p, ...
               'pAIC',pAIC,'aic',aic,'shares',shares,'MSPE',MSPE,'D',D, ...
               'table1',chk,'seasonality',seas,'Xexo',Xexo);

    if opts.gate && ~chk.pass
        fprintf(2, ['Stopping before the scenarios: the estimated VAR does not ' ...
                    'reproduce Table 1.\n' ...
                    'The returned struct therefore has NO .Fscen, .scen, .Fbase, ' ...
                    '.Fbk or .wzinfo\nfields -- only the model, the decomposition ' ...
                    'and the Table 1 check. Re-run with\nopts.gate = false to ' ...
                    'compute the scenarios and draw the figures anyway.\n\n']);
        return;
    end

    % ================= 4. baseline and WZ conditional forecasts ==========
    Ylast = Y(end-p+1:end, :);
    Fbase = baselineForecast(model, Ylast, H);

    scen   = buildScenarios(tr, H, opts.scenarios);
    nS     = numel(scen);
    Fscen  = cell(1,nS);  wzinfo = cell(1,nS);
    fprintf('Waggoner-Zha conditional forecasts:\n');
    for k = 1:nS
        sdk = [];
        if isfield(scen(k),'sd'), sdk = scen(k).sd; end
        [Fc, ~, info] = wzConditional(model, Fbase, scen(k).cvar, ...
                                      scen(k).path, H, MAd, [], sdk);
        Fscen{k} = Fc;  wzinfo{k} = info;
        % ||eps|| against its expectation under the model: eps has identity
        % covariance, so a typical draw has E||eps||^2 = K*H. A scenario needing
        % much more than that is far out in the model's tail, which is worth
        % knowing before reading magnitudes off the figure.
        fprintf('  %-24s residual = %.2e   ||eps|| = %6.2f  (%.1fx typical)\n', ...
                scen(k).name, info.constraintResidual, info.normEps, ...
                info.normEps / sqrt(model.K*H));
    end

    % ---- BK robustness for the monetary-policy scenario -----------------
    % NOTE: this is NOT the paper's Figure 4 exercise. The paper re-estimates
    % on a pre-1992 sample, picks the 48 months of 1988-1991 with the largest
    % policy-shock variability, and feeds those structural shocks to BK. That
    % needs data this package does not load. What follows takes the WZ
    % innovations, maps them to structural space, and keeps only the policy
    % shock -- the "tightening driven purely by policy shocks" reading.
    % The BK robustness check needs a policy-rate scenario; find one rather
    % than assuming it is the fourth, since the spec may be custom.
    kMP = [];
    for k = 1:nS
        if any(scen(k).cvar == idx.rate), kMP = k; break; end
    end
    if isempty(kMP), kMP = 1; end
    Emp  = zeros(model.K, H);
    Emp(idx.rate, :) = wzinfo{kMP}.eps(idx.rate, :);
    Fbk  = bkConditional(model, Fbase, Emp, H, D, MAd);
    fprintf(['  BK keeps the policy shock only: ||eps_mp|| = %.2f of ' ...
             '||eps|| = %.2f\n'], norm(Emp(:)), wzinfo{kMP}.normEps);

    R.Fbase = Fbase;  R.Fscen = Fscen;  R.scen = scen;  R.Fbk = Fbk;
    R.wzinfo = wzinfo;

    % ================= uncertainty bands (optional) ======================
    if opts.bands
        fprintf('\nUncertainty bands (shock + parameter)...\n');
        bo = opts.bandOpts;
        if ~isempty(priorStruct) && ~isfield(bo,'prior'), bo.prior = priorStruct; end
        R.bands = wzBands(Y, p, Xexo, scen, H, bo);
        fprintf(['  %d draws (%s), %d explosive replications redrawn.\n' ...
                 '  Read scen(k).diff for whether a scenario is ' ...
                 'distinguishable from the baseline:\n' ...
                 '  it carries parameter uncertainty only, since the shock ' ...
                 'uncertainty common\n  to the scenario and the baseline ' ...
                 'cancels in the difference.\n'], ...
                R.bands.nDraws, R.bands.bandSource, R.bands.nRejected);
    end

    % ================= figures ===========================================
    if ~opts.figures, return; end
    fprintf('\nRendering figures to %s ...\n', outdir);
    try
        xlo = datenum(opts.figStart(1), opts.figStart(2), 1);
        dvL = datevec(data.dates(end));
        nL  = dvL(1)*12 + (dvL(2)-1) + H;
        xhiMax = datenum(floor(nL/12), mod(nL,12)+1, 1);
        if isempty(opts.figEnd)
            xhi = xhiMax;
        else
            xhi = datenum(opts.figEnd(1), opts.figEnd(2), 1);
            if xhi > xhiMax
                fprintf(2, ['figEnd %s is past the forecast horizon (%s); ' ...
                            'clipped.\n'], datestr(xhi,'yyyy-mm'), ...
                            datestr(xhiMax,'yyyy-mm'));
                xhi = xhiMax;
            elseif xhi <= data.dates(end)
                error('main_risk_scenarios:figEnd', ...
                    'figEnd %s is at or before the forecast origin %s.', ...
                    datestr(xhi,'yyyy-mm'), datestr(data.dates(end),'yyyy-mm'));
            end
        end
        nMade = 0;
        if opts.paperFigs
            figure1_scenarios(data, idx, Fbase, Fscen, scen, H, outdir, xlo, xhi);
            figure2_canada(data, idx, Fbase, Fscen, scen, H, outdir, xlo, xhi);
            figure3_levels(data, idx, Fbase, Fscen, scen, H, outdir, xlo, xhi);
            figure4_mp_robust(data, idx, Fbase, Fscen{kMP}, Fbk, H, outdir, xlo, xhi);
            nMade = nMade + 4;
        end
        if opts.bands && isfield(R,'bands') && isfield(R.bands,'condDraws')
            for k = 1:nS
                figureBands(data, idx, Fbase, R.bands, k, H, outdir, xlo, [], xhi);
                nMade = nMade + 1;
                if opts.diffFigs
                    figureDiffBands(data, idx, R.bands, k, H, outdir, xlo, xhi);
                    nMade = nMade + 1;
                end
            end
        end
        if nMade == 0
            fprintf(2, ['No figures requested: paperFigs is false and either ' ...
                        'bands is false or\nthe draws are missing.\n']);
        else
            fprintf('Figures saved (%d).\n', nMade);
        end
    catch ME
        fprintf(2,'Plotting skipped (%s). Numeric results are still returned.\n', ME.message);
    end
end

% =====================================================================
function o = withDefaults(o, d)
    f = fieldnames(d);
    for i = 1:numel(f)
        if ~isfield(o, f{i}), o.(f{i}) = d.(f{i}); end
    end
end

% =====================================================================
%                              FIGURES
% =====================================================================
% Display conventions follow the paper's axes: every growth-rate panel is a
% 12-month percent change. Annualised monthly growth (x1200) was used here
% previously and is far too volatile to match -- April 2020 U.S. industrial
% production alone annualises to about -160%, while the paper's panel runs
% over +/-20.
% =====================================================================

function figure1_scenarios(data, idx, Fbase, Fscen, scen, H, outdir, xlo, xhi)
% Figure 1 analogue: each panel shows the targeted variable under history
% (green), baseline (black) and its own scenario (colour).
    f = figure('Visible','off','Position',[100 100 900 650]);
    specs = { idx.oil,  'level', 'High Oil Prices',        1;
              idx.usip, 'yoy',   'U.S. Recession',         2;
              idx.unemp,'level', 'Low Unemployment',       3;
              idx.rate, 'level', 'Tight Monetary Policy',  4};
    for i = 1:4
        subplot(2,2,i);
        k = specs{i,4};
        drawVar(data, specs{i,1}, specs{i,2}, Fbase, {Fscen{k}}, ...
                {scen(k).color}, H, specs{i,3}, xlo, xhi);
    end
    saveFig(f, fullfile(outdir,'figure1_scenarios.png'));
end

function figure2_canada(data, idx, Fbase, Fscen, scen, H, outdir, xlo, xhi)
% Figure 2 analogue: Canadian variables, baseline + all four scenarios.
%
% Colours follow the paper's Section 5 TEXT: purple = oil, orange = U.S.
% recession, blue = unemployment, yellow = monetary policy. Figure 3's legend
% in the paper agrees; Figure 2's published legend rotates three of the four
% labels and should be ignored when comparing curve by curve.
    f = figure('Visible','off','Position',[100 100 1000 750]);
    specs = { idx.cpi,  'yoy',   'Prices (YoY %)';
              idx.house,'yoy',   'Housing Prices (YoY %)';
              idx.unemp,'level', 'Unemployment (%)';
              idx.rate, 'level', 'BC Rate (%)';
              idx.fx,   'yoy',   'Exchange Rate (YoY %)'};
    for i = 1:5
        subplot(3,2,i);
        drawVar(data, specs{i,1}, specs{i,2}, Fbase, Fscen, {scen.color}, H, specs{i,3}, xlo, xhi);
    end
    subplot(3,2,6); axis off; hold on;
    plot(nan,nan,'k-','LineWidth',1.8);
    for k=1:numel(scen), plot(nan,nan,'-','Color',scen(k).color,'LineWidth',1.6); end
    legend([{'Baseline'}, {scen.name}], 'Location','west','Box','off');
    saveFig(f, fullfile(outdir,'figure2_canada.png'));
end

function figure3_levels(data, idx, Fbase, Fscen, scen, H, outdir, xlo, xhi)
% Figure 3 analogue: housing and FX in levels.
    f = figure('Visible','off','Position',[100 100 950 380]);
    subplot(1,2,1);
    drawVar(data, idx.house,'level', Fbase, Fscen, {scen.color}, H, 'Housing Prices (index)', xlo, xhi);
    subplot(1,2,2);
    drawVar(data, idx.fx,'level', Fbase, Fscen, {scen.color}, H, 'Exchange Rate (CAD/USD)', xlo, xhi);
    legend([{'Baseline'}, {scen.name}], 'Location','best','Box','off');
    saveFig(f, fullfile(outdir,'figure3_levels.png'));
end

function figure4_mp_robust(data, idx, Fbase, Fwz, Fbk, H, outdir, xlo, xhi)
% Figure 4 analogue: monetary-policy tightening, WZ vs BK.
    f = figure('Visible','off','Position',[100 100 1000 750]);
    specs = { idx.cpi,  'yoy',   'Prices (YoY %)';
              idx.house,'yoy',   'Housing Prices (YoY %)';
              idx.unemp,'level', 'Unemployment (%)';
              idx.rate, 'level', 'BC Rate (%)';
              idx.fx,   'yoy',   'Exchange Rate (YoY %)'};
    cWZ = [0.93 0.69 0.13]; cBK = [0.49 0.42 0.75];
    for i = 1:5
        subplot(3,2,i);
        drawVar(data, specs{i,1}, specs{i,2}, Fbase, {Fwz, Fbk}, {cWZ, cBK}, H, specs{i,3}, xlo, xhi);
    end
    subplot(3,2,6); axis off; hold on;
    plot(nan,nan,'k-','LineWidth',1.8); plot(nan,nan,'-','Color',cWZ,'LineWidth',1.6);
    plot(nan,nan,'-','Color',cBK,'LineWidth',1.6);
    legend({'Baseline','Waggoner-Zha','Baumeister-Kilian'},'Location','west','Box','off');
    saveFig(f, fullfile(outdir,'figure4_mp_robustness.png'));
end

% =====================================================================
%                          PLOT PRIMITIVES
% =====================================================================
function drawVar(data, v, mode, Fbase, Fscens, colors, H, ttl, xlo, xhi)
% xlo (optional) is the left-hand x limit as a datenum; default is 48 months
% of history, i.e. the same span the paper shows.
    last = data.dates(end);
    [bx, bval] = dispSeries(data, v, mode, Fbase(:,v), H);
    hist = bx <= last;  fut = bx >= last;
    hold on; box on;
    plot(bx(hist), bval(hist), '-', 'Color',[0.15 0.6 0.25], 'LineWidth',1.3);
    plot(bx(fut),  bval(fut),  'k-', 'LineWidth',1.7);
    for k = 1:numel(Fscens)
        [sx, sval] = dispSeries(data, v, mode, Fscens{k}(:,v), H);
        sfut = sx >= last;
        plot(sx(sfut), sval(sfut), '-', 'Color',colors{k}, 'LineWidth',1.4);
    end
    yl = ylim; plot([last last], yl, '-', 'Color',[0.6 0.6 0.6]);
    title(ttl); grid on;
    if nargin < 9  || isempty(xlo), xlo = addMonths(last,-47); end
    if nargin < 10 || isempty(xhi), xhi = addMonths(last,H);   end
    xlim([xlo, xhi]);
    datetick('x','yyyy','keeplimits');
end

function dn = addMonths(d, k)
% Month arithmetic done on a month counter rather than by handing datenum an
% out-of-range month. datenum(2022, 12-47, 1) does NOT reliably roll back into
% 2019 -- it lands on 2022-01 -- which silently cropped the plots to twelve
% months of history instead of forty-eight.
    dv = datevec(d);
    n  = dv(1)*12 + (dv(2)-1) + k;
    dn = datenum(floor(n/12), mod(n,12)+1, 1);
end

function saveFig(f, path)
    try
        print(f, path, '-dpng', '-r120');
    catch
        saveas(f, path);
    end
    close(f);
end
