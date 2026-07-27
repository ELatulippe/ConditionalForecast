function outFile = build_dashboard(opts)
%BUILD_DASHBOARD  Build the self-contained interactive dashboard (dashboard.html).
%
%   build_dashboard
%   build_dashboard(struct('nDraws',400,'open',true))
%
% Estimates the VAR (OFFLINE, from cache/panel_k9.mat), computes the structural
% impulse responses and their bootstrap/posterior bands, embeds them together
% with the scenario fan-chart PNGs already in figures/, and writes a single
% self-contained file:
%
%     dashboard/dashboard.html
%
% Open that file in any browser (no server, no internet). Three tabs:
%   * Impulse-response explorer -- recomputes the IRFs live as you change the
%     shock, its size, the horizon and the cumulate/band toggles.
%   * Scenario builder (live)   -- build a Waggoner-Zha conditional forecast in
%     the browser (ramp/hold/pulse/policy/replay constraints, joint & soft), and
%     see the baseline-vs-scenario level paths recompute as you drag the inputs.
%   * Scenario gallery          -- flip through the precomputed scenario charts.
%
% OPTS (all optional)
%   .nDraws  bootstrap/posterior draws for the IRF bands (default 400; 0 = none)
%   .config  options struct for the panel/estimate (default scenario_config())
%   .open    try to open the file in the browser when done (default false)
%
% Requires MATLAB R2016b+ (uses jsonencode / matlab.net.base64encode). Run the
% scenario scripts (or run_all) first if you want the gallery populated.

    if nargin < 1, opts = struct(); end
    if ~isfield(opts,'nDraws'), opts.nDraws = 400; end
    if ~isfield(opts,'open'),   opts.open   = false; end

    root = fileparts(mfilename('fullpath'));
    addpath(root); addpath(fullfile(root,'src'));
    assert(exist('main_risk_scenarios','file')==2, ...
        'main_risk_scenarios.m not found under ./src -- keep this script at the project root.');
    assert(exist('jsonencode','file')>0 || exist('jsonencode','builtin')>0, ...
        'build_dashboard needs MATLAB R2016b+ (jsonencode).');

    tplFile = fullfile(root,'dashboard','dashboard_template.html');
    assert(exist(tplFile,'file')==2, 'Missing template: %s', tplFile);

    % ---- 1. estimate offline -------------------------------------------
    if isfield(opts,'config') && ~isempty(opts.config), ob = opts.config;
    else,                                               ob = scenario_config();
    end
    ob.offline = true;  ob.bands = false;  ob.figures = false;  ob.gate = false;
    fprintf('build_dashboard: estimating (offline)...\n');
    R = main_risk_scenarios('fred', ob);

    names = R.data.names;  tcode = R.data.tcode;  K = R.model.K;  H = ob.H;

    % ---- 2. impulse responses: raw AND cumulated, each with bands ------
    fprintf('build_dashboard: impulse responses (%d draws)...\n', opts.nDraws);
    optRaw = struct('cumulate',false,'scaleTo',[],'boot',opts.nDraws,'verbose',false);
    optCum = struct('cumulate',true, 'scaleTo',[],'boot',opts.nDraws,'verbose',false);
    Iraw = irfCholesky(R, H, optRaw);
    Icum = irfCholesky(R, H, optCum);
    hasBand = isfield(Iraw,'lo') && isfield(Icum,'lo');

    data = struct();
    [dScale, dUnit] = displayScales(names, tcode, R.data.levels(end,:));
    data.meta = struct('K',K, 'H',H, 'names',{names}, 'tcode',{tcode}, ...
                       'key',{keyVars(names)}, ...
                       'dispScale',dScale, 'dispUnit',{dUnit}, ...
                       'generated', datestr(now,'yyyy-mm-dd HH:MM'), ...
                       'subtitle', 'Conditional-forecast (Waggoner-Zha) replication');
    data.irf = struct('hasBand', hasBand, ...
                      'band_q', Iraw.quantiles, ...
                      'raw', modeStruct(Iraw, hasBand), ...
                      'cum', modeStruct(Icum, hasBand));

    % ---- 2b. forecast model for the LIVE scenario builder --------------
    % Fbase (baseline forecast, transformed units), the last observed levels,
    % and the full raw-level history + dates (for 'replay' and the chart tail).
    p     = R.model.p;
    Ylast = R.tr.Y(end-p+1:end, :);
    Fbase = baselineForecast(R.model, Ylast, H);
    data.model = struct( ...
        'Fbase',      round(Fbase*1e6)/1e6, ...                 % H x K -> [h][i]
        'lastLevels', round(R.data.levels(end,:)*1e4)/1e4, ...  % 1 x K
        'levels',     round(R.data.levels*1e4)/1e4, ...         % T x K -> [t][i]
        'dates',      {cellstr(datestr(R.data.dates,'yyyy-mm'))}, ...
        'A',          round(R.model.A*1e6)/1e6, ...             % Kp x Kp companion
        'mu',         round(R.model.mu(:).'*1e6)/1e6, ...       % 1 x Kp
        'p',          R.model.p, ...
        'Y',          round(R.tr.Y*1e6)/1e6, ...                % (T-1) x K transformed history
        'trDates',    {cellstr(datestr(R.tr.dates,'yyyy-mm'))});

    % ---- 3. scenario gallery from existing figures ---------------------
    data.scenarios = gatherScenarios(fullfile(root,'figures'));

    % ---- 4. encode + inject into the template --------------------------
    js = jsonencode(data);
    tpl = fileread(tplFile);
    marker = '/*__DATA__*/ null';
    if isempty(strfind(tpl, marker))                                   %#ok<STREMP>
        error('build_dashboard:tpl','Template marker not found in %s', tplFile);
    end
    html = strrep(tpl, marker, ['/*__DATA__*/ ' js]);

    outFile = fullfile(root,'dashboard','dashboard.html');
    fid = fopen(outFile,'w');
    assert(fid>=0, 'Cannot write %s', outFile);
    fwrite(fid, html);  fclose(fid);

    fprintf('build_dashboard: wrote %s (%.1f MB, %d scenarios).\n', ...
            outFile, numel(html)/1e6, numel(data.scenarios));
    fprintf('Open it in a browser: the IRF tab recomputes live, no server needed.\n');

    if opts.open
        try
            if ispc, winopen(outFile);
            elseif ismac, system(['open "' outFile '"']);
            else, system(['xdg-open "' outFile '" &']);
            end
        catch, end
    end
end

% ======================================================================
function m = modeStruct(I, hasBand)
% Pack one IRF mode (raw or cumulated) into nested [i][j][h] cells.
    m = struct('point', {pack3(I.resp)});
    if hasBand
        m.lo = pack3(I.lo);
        m.hi = pack3(I.hi);
    end
end

function C = pack3(A)
% A is K x K x H  ->  C{i}{j} = 1 x H row (rounded to keep the file small).
    [K1,K2,Hh] = size(A);
    C = cell(K1,1);
    for i = 1:K1
        C{i} = cell(1,K2);
        for j = 1:K2
            v = reshape(A(i,j,:), 1, Hh);
            C{i}{j} = round(v*1e5)/1e5;      % round(x,5) without newer syntax
        end
    end
end

function [scale, unit] = displayScales(names, tcode, lastLev)
% Per-variable DISPLAY scaling for the level charts (does not touch the model).
% Employment (FRED LFEMTTTTCAM647S) is in persons; show it in thousands. The
% chosen scale adapts to the magnitude so it reads cleanly whatever the vintage.
    K = numel(names);
    scale = ones(1,K);
    unit  = repmat({''}, 1, K);
    for j = 1:K
        if strcmpi(names{j}, 'emp')
            lv = abs(lastLev(j));
            if     lv >= 1e6, scale(j) = 1e-3; unit{j} = 'level, 000 persons';   % persons -> thousands
            elseif lv >= 1e3, scale(j) = 1;    unit{j} = 'level, 000 persons';   % already thousands
            else,             scale(j) = 1;    unit{j} = 'level';
            end
        end
    end
end

function key = keyVars(names)
% A small default set to show first, intersected with what the panel has.
    want = {'oil','cpi','unemp','rate','fx','house'};
    key = want(ismember(want, names));
    if isempty(key), key = names(1:min(6,numel(names))); end
end

% ======================================================================
function S = gatherScenarios(figdir)
% Collect scenario fan-chart PNGs (embedded as data URIs) for the gallery,
% GROUPED by kind, each with a curated title and a one-line description.
    secs = { ''                , 'Single-variable scenarios'; ...
             'joint_scenarios' , 'Joint (multi-variable) scenarios'; ...
             'temporary_shocks', 'Temporary / partial-horizon shocks'; ...
             'soft_conditions' , 'Soft (distributional) conditions'; ...
             'counterfactual'  , 'Counterfactuals' };
    S = struct('title',{}, 'group',{}, 'desc',{}, 'img',{});
    for r = 1:size(secs,1)
        sub = secs{r,1};  grp = secs{r,2};
        d = figdir; if ~isempty(sub), d = fullfile(figdir, sub); end
        if ~exist(d,'dir'), continue; end
        L = dir(fullfile(d,'*.png'));
        L = L(~[L.isdir]);
        [~,ord] = sort({L.name});
        for k = ord(:)'
            p = fullfile(d, L(k).name);
            [ttl, dsc] = scenarioMeta(L(k).name);
            S(end+1) = struct('title',ttl, 'group',grp, 'desc',dsc, ...
                              'img',imgURI(p,figdir)); %#ok<AGROW>
        end
    end
end

function [ttl, dsc] = scenarioMeta(fname)
% Curated title + one-line description, matched by a SUBSTRING of the file name
% (robust to the run-time numbers baked into some names, e.g. oil_30_to_133).
% Order matters: more specific keys before their prefixes.
    key = lower(fname);
    rules = { ...
      'oil_30_to',                'Oil +30%',                          'WTI ramps ~30% above the forecast-origin price over three months, then holds.'; ...
      'tight_labour_market',      'Tight labour market, 1yr',          'Unemployment held ~1.5pp below its latest value for a year, then free.'; ...
      'gradual_tightening',       'Gradual tightening (+150bp)',       'Policy rate stepped up 25bp at a time to +150bp, held a year, then eased.'; ...
      'u_s_recession_boc_eases',  'U.S. recession, BoC eases',         'US industrial production replays 2008-09 while the Bank cuts the policy rate ~200bp.'; ...
      'u_s_recession_no_policy',  'U.S. recession, no policy response','US activity down 2008-style with the policy rate held fixed -- a policy counterfactual.'; ...
      'u_s_recession',            'U.S. recession',                    'US industrial production replays its 2008-09 recession path.'; ...
      'oil_boom',                 'Oil boom + tightening',             'Oil rises ~40% and the Bank leans against it with gradual tightening.'; ...
      'risk_off',                 'Risk-off',                          'Equity crash (2008 replay), flight-to-quality in yields and a weaker CAD, together.'; ...
      'u_s_reflation',            'U.S. reflation',                    'US CPI replays its 2021 surge and the T-bill path rises ~200bp; Canada responds endogenously.'; ...
      'stagflation',              'Stagflation',                       'Oil up and unemployment up together -- forcing the system onto a supply-shock path.'; ...
      'oil_spike_loosely',        'Oil spike, loosely held',           'An oil pulse conditioned as a distribution rather than an exact path.'; ...
      'oil_spike_6',              'Oil spike, 6 months',               'Oil pulses to $110 over three months, holds three, then the model sets the tail.'; ...
      'oil_at_110',               'Oil at $110, 6m then free',         'Oil pinned at $110 for six months only; the level then evolves freely.'; ...
      'rate_150bp',               'Rate +150bp, 1yr then free',        'A temporary +150bp tightening for a year that then normalises endogenously.'; ...
      'rate_100bp',               'Rate +100bp, tight then loose',     '+100bp with a per-horizon sd -- tightly held early, loosely later.'; ...
      'tight_labour_9m',          'Tight labour, 9m then free',        'Unemployment held low for nine months only, then free.'; ...
      'cad_weakens',              'CAD weakens, 1Q then free',         'A one-quarter CAD depreciation to ~1.45, then free.'; ...
      'equity_selloff',           'Equity selloff, 6m then free',      'A six-month replay of a 2008 equity drawdown, then recovery left to the model.'; ...
      'tightening_roughly',       'Tightening, roughly',               'A tightening path conditioned softly (sd 0.4pp) so the rate itself keeps a band.'; ...
      'unemp_5',                  'Unemp ~5%, 1yr (soft)',             'Unemployment around 5% for a year, softly (sd 0.5pp).'; ...
      'cad_1_42',                 'CAD ~1.42 (soft)',                  'CAD around 1.42, softly conditioned.'; ...
      'tighter_boc',              'Tighter BoC (+100bp)',              'What if the Bank had held the policy rate 100bp higher over 2022-23? Price level and unemployment.'; ...
      'no_oil_spike',             'No 2022 commodity spike',           'What if the 2022 oil spike had not happened? Price level and unemployment.' };
    for i = 1:size(rules,1)
        if ~isempty(strfind(key, rules{i,1}))            %#ok<STREMP>
            ttl = rules{i,2};  dsc = rules{i,3};  return;
        end
    end
    ttl = humanize(fname);  dsc = '';
end

function uri = imgURI(pngPath, figdir)
% Base64 data URI so the dashboard is self-contained; relative path fallback.
    try
        fid = fopen(pngPath,'r');  raw = fread(fid, Inf, '*uint8');  fclose(fid);
        uri = ['data:image/png;base64,' char(matlab.net.base64encode(raw))];
    catch
        rel = strrep(strrep(pngPath, figdir, '../figures'), '\', '/');
        uri = rel;
    end
end

function s = humanize(name)
    [~, name] = fileparts(name);
    for pre = {'figure5_bands_','figure_','fevd_','irf_'}
        if strncmp(name, pre{1}, numel(pre{1})), name = name(numel(pre{1})+1:end); end
    end
    s = strtrim(strrep(name, '_', ' '));
    % fix common acronyms when they appear as whole words
    s = regexprep(s, '\<u s\>',  'US');
    s = regexprep(s, '\<us\>',   'US');
    s = regexprep(s, '\<cad\>',  'CAD');
    s = regexprep(s, '\<boc\>',  'BoC');
    s = regexprep(s, '\<goc\>',  'GoC');
    s = regexprep(s, '\<tsx\>',  'TSX');
    s = regexprep(s, '\<cpi\>',  'CPI');
    if ~isempty(s), s(1) = upper(s(1)); end
end
