function [data, span] = loadCanadaData(startY, startM, endY, endM, source, opts)
% LOADCANADADATA  Assemble the 7-variable monthly panel used in the paper.
%
%   data = loadCanadaData(1992,1, 2022,12, 'fred', opts)
%   data = loadCanadaData(1992,1, 2022,12, 'local', opts)
%   data = loadCanadaData(1992,1, 2022,12, 'synthetic')
%
% Variable order (as they enter y_t; ordering matters for Cholesky):
%   1 oil    WTI crude, USD/bbl              MCOILWTICO      -> dlog
%   2 usip   US industrial production        INDPRO          -> dlog
%   3 cpi    Canada CPI, all items           CPALCY01CAM661N -> dlog
%   4 house  Canada house-price index        user-supplied   -> dlog
%   5 unemp  Canada unemployment rate, %     LRUNTTTTCAM156S -> level
%   6 rate   Canada 3-month interbank, %     IR3TIB01CAM156N -> level
%   7 fx     CAD per USD                     EXCAUS          -> dlog
%
% OPTS (struct, all fields optional)
%   .houseFile   path to a monthly CSV holding the house-price index.
%                REQUIRED for 'fred' and 'local' unless .allowHouseSpline is
%                set -- see the housing note below. Either layout works:
%                a long two-column 'date,value' file, or Statistics Canada's
%                "Download as displayed" export with months across a header
%                row (see readSeriesCSV).
%   .houseSeries which row to take from a wide StatCan export. Default
%                'Total (house and land)'. Use e.g. 'House only' to take the
%                structure-only component.
%   .allowHouseSpline  logical, default false. If true, fall back to the
%                quarterly BIS series QCAN628BIS spline-interpolated to
%                monthly. This is a diagnostic escape hatch only: it makes
%                monthly housing growth ~95% predictable from its own lags,
%                collapses that equation's innovation variance, and thereby
%                distorts Sigma_u, the Cholesky factor D, every Theta_j and
%                all four conditional forecasts. Do not use it to replicate.
%   .csvDir      directory of pre-downloaded FRED CSVs (source = 'local')
%   .cacheFile   .mat file to save/load the assembled panel (source='fred')
%   .cpiFile     path to a seasonally adjusted CPI CSV, overriding the FRED
%                default CPALCY01CAM661N, which is NOT seasonally adjusted.
%                Statistics Canada table 18-10-0006-01 (Consumer Price Index,
%                monthly, seasonally adjusted) downloads in the same wide
%                layout readSeriesCSV handles.
%   .cpiSeries   which row to take from a wide CPI export. Default tries
%                'All-items' / 'Ensemble'.
%   .rateFile    path to a monthly policy/short-rate CSV, overriding the FRED
%                default IRSTCB01CAM156N (the Bank Rate). The Bank Rate is a
%                step function that changes in discrete 25bp jumps and behaves
%                close to exogenously; a market rate such as the 3-month
%                treasury bill responds far more to macro developments, which
%                matters for how much of the rate's forecast-error variance is
%                its own shock at long horizons.
%   .rateSeries  which row to take from a wide rate export.
%   .files       struct keyed by variable name ('oil','usip','cpi','house',
%                'unemp','rate','fx') giving a CSV path for any series you
%                want to supply yourself. This supersedes houseFile/cpiFile/
%                rateFile, which still work and map into it. Use it when a
%                FRED series stops updating and caps the panel:
%                   opts.files = struct('unemp','1410028701-eng.csv');
%   .rows        struct keyed the same way, giving the row label to take from
%                a wide "download as displayed" export.
%   .extra       cell array of structs, each ADDING a variable to the panel:
%                  .name   short name, becomes a field of data.idx
%                  .id     FRED series id,  OR
%                  .file   CSV path (+ .row for a wide export)
%                  .tcode  'dlog' or 'level'  (default 'level')
%                  .after  name of the variable it should FOLLOW in the
%                          recursive ordering (default: last)
%                Ordering is identification, not bookkeeping: a variable
%                placed at position k is assumed not to respond within the
%                month to anything after it. Financial variables that price
%                continuously usually belong late, near the exchange rate.
%                Each extra variable costs K*p+1 more coefficients per
%                equation and adds an equation, so past nine or so variables
%                OLS on ~410 months gets thin and a shrinkage prior is the
%                honest next step.
%   .ids         struct keyed by variable name overriding the FRED series ID,
%                e.g. struct('unemp','LRUN64TTCAM156S'). Lower friction than
%                .files when the replacement also lives on FRED: nothing to
%                download. The OECD-sourced defaults for unemployment and the
%                policy rate are the ones that go stale.
%
% For backward compatibility a char 6th argument is read as .cacheFile for
% source 'fred' and as .csvDir for source 'local'.
%
% HOUSING NOTE
%   The paper's real-estate series is the monthly Canadian house-price index
%   in the Fortin-Gagnon et al. (2022) database, whose 12-month growth is
%   about 0% in 2019, peaks near 11.5% in late 2021, and is about 4% at the
%   end of 2022 (compare Figure 2). That profile matches Statistics Canada's
%   New Housing Price Index (table 18-10-0205), not the BIS residential
%   property price index, whose 12-month growth peaks around 26% over the
%   same window. FRED carries no monthly Canadian house-price index, so
%   export the monthly index to a two-column CSV and pass it via
%   opts.houseFile.
%
% OUTPUT
%   data  struct with .dates (T x 1 datenum), .levels (T x 7), .names,
%         .tcode, .idx
%   span  struct array, one entry per series, with .name, .first, .last
%         (datenums of the first/last observed month) and .binding, a flag
%         marking the series that truncated the common sample.

    if nargin < 5 || isempty(source), source = 'fred'; end
    if nargin < 6, opts = struct(); end
    if ischar(opts)                      % legacy positional 6th argument
        if strcmpi(source,'local'), opts = struct('csvDir', opts);
        else,                       opts = struct('cacheFile', opts);
        end
    end
    opts = withDefaults(opts, struct('houseFile','', 'houseSeries','', ...
                                     'allowHouseSpline',false, ...
                                     'csvDir','', 'cacheFile','', ...
                                     'dataDir','', 'cacheDir','', ...
                                     'offline',false, ...
                                     'cpiFile','', 'cpiSeries','', ...
                                     'rateFile','', 'rateSeries','', ...
                                     'files',struct(), 'rows',struct(), ...
                                     'ids',struct(), 'extra',{{}}, ...
                                     'openEnded',false));

    names = {'oil','usip','cpi','house','unemp','rate','fx'};
    tcode = {'dlog','dlog','dlog','dlog','level','level','dlog'};
    % IR3TIB01CAM156N (3-month interbank) replaces IRSTCB01CAM156N (the Bank
    % Rate), which stopped updating: FRED's migration to the new OECD Data
    % Explorer carried only some Main Economic Indicators series across. The
    % test is the Notes block on the series page -- "OECD Data Filters:
    % REF_AREA: ..." means migrated and live, the older "OECD Descriptor ID:"
    % style means frozen. This is a CONCEPT CHANGE, not just a longer series:
    % the Bank Rate is a policy instrument that moves in 25bp steps, while the
    % interbank rate is market-determined and embeds expectations of policy
    % over the coming quarter. Revert with opts.ids.rate = 'IRSTCB01CAM156N'.
    ids   = {'MCOILWTICO','INDPRO','CPALCY01CAM661N','QCAN628BIS', ...
             'LRUNTTTTCAM156S','IR3TIB01CAM156N','EXCAUS'};
    quarterly = [false false false true false false false];

    % ---- splice in any extra variables ---------------------------------
    extra = opts.extra;
    if isstruct(extra), extra = num2cell(extra); end
    extraFile = repmat({''}, 1, numel(names));
    extraRow  = repmat({''}, 1, numel(names));
    for e = 1:numel(extra)
        ex = extra{e};
        if ~isfield(ex,'name') || isempty(ex.name)
            error('loadCanadaData:extra', 'Each .extra entry needs a .name.');
        end
        if any(strcmp(ex.name, names))
            error('loadCanadaData:extra', ...
                'Variable "%s" is already in the panel.', ex.name);
        end
        if isfield(ex,'after') && ~isempty(ex.after)
            pos = find(strcmp(ex.after, names), 1);
            if isempty(pos)
                error('loadCanadaData:extra', ...
                    '.after = "%s" is not a variable in the panel (%s).', ...
                    ex.after, strjoin(names, ', '));
            end
        else
            pos = numel(names);
        end
        exId   = ''; if isfield(ex,'id'),   exId   = ex.id;   end
        exFile = ''; if isfield(ex,'file'), exFile = ex.file; end
        exRow  = ''; if isfield(ex,'row'),  exRow  = ex.row;  end
        exT    = 'level'; if isfield(ex,'tcode') && ~isempty(ex.tcode), exT = ex.tcode; end
        if isempty(exId) && isempty(exFile)
            error('loadCanadaData:extra', ...
                'Variable "%s" needs either .id (FRED) or .file (CSV).', ex.name);
        end
        ins = @(c, v) [c(1:pos), {v}, c(pos+1:end)];
        names     = ins(names, ex.name);
        tcode     = ins(tcode, exT);
        ids       = ins(ids,   exId);
        extraFile = ins(extraFile, exFile);
        extraRow  = ins(extraRow,  exRow);
        quarterly = [quarterly(1:pos), false, quarterly(pos+1:end)];
    end

    nV = numel(names);
    idx = struct();
    for j = 1:nV, idx.(names{j}) = j; end

    % Per-variable FRED ID overrides.
    for j = 1:nV
        if isfield(opts.ids, names{j}) && ~isempty(opts.ids.(names{j}))
            ids{j} = opts.ids.(names{j});
            quarterly(j) = false;          % assume monthly unless told otherwise
            fprintf('  %-6s FRED id overridden -> %s\n', names{j}, ids{j});
        end
    end

    grid = monthGrid(startY,startM,endY,endM);
    T = numel(grid);  L = nan(T, nV);

    switch lower(source)
    case 'synthetic'
        L = syntheticPanel(T);
        if nV > 7                      % fill any extra columns with noise
            randn('seed',7);
            L = [L, 100*exp(0.001*(1:T).' + 0.02*cumsum(randn(T, nV-7)))];
        end
        [data, span] = pack(grid, L, names, tcode, idx, false, false);
        return;

    case {'fred','local'}
        cachePath = resolveCache(opts.cacheDir, opts.cacheFile);
        if strcmpi(source,'fred') && ~isempty(cachePath) && exist(cachePath,'file')
            S = load(cachePath);
            % The cache is read before anything else, so a stale one would
            % silently override .extra and hand back a panel with the wrong
            % variables. Check the variable set before trusting it.
            if isequal(S.data.names(:).', names(:).')
                data = S.data;  span = seriesSpans(data.dates, data.levels, names);
                return;
            end
            if opts.offline
                error('loadCanadaData:offline', ...
                    ['Offline mode, but the cache %s holds a DIFFERENT variable set.\n' ...
                     '    cached: %s\n    wanted: %s\n' ...
                     'Rebuild the cache for this variable set with ONE online run\n' ...
                     '(delete the file and run freeze_panel, or main_risk_scenarios\n' ...
                     'with source=''fred'' and internet), then re-run offline.'], ...
                    cachePath, strjoin(S.data.names, ' '), strjoin(names, ' '));
            end
            fprintf(2, ['\n*** Cache %s holds a different variable set and is ' ...
                        'being IGNORED.\n    cached: %s\n    wanted: %s\n' ...
                        '    Refetching; the file will be overwritten. ***\n\n'], ...
                    opts.cacheFile, strjoin(S.data.names, ' '), strjoin(names, ' '));
        end
        % Offline mode never touches the network: a usable cache must exist.
        if opts.offline && strcmpi(source,'fred')
            error('loadCanadaData:offline', ...
                ['Offline mode: no usable panel cache at\n    %s\n' ...
                 'Nothing can be fetched from FRED. Do ONE online run to build it\n' ...
                 '(run freeze_panel, or call main_risk_scenarios with source=''fred''\n' ...
                 'and internet access), then re-run offline. For a run entirely from\n' ...
                 'pre-downloaded FRED CSVs instead, use source=''local'' with a csvDir.'], ...
                cachePath);
        end
        haveHouse = ~isempty(opts.houseFile) || ...
                    (isfield(opts.files,'house') && ~isempty(opts.files.house));
        if ~haveHouse && ~opts.allowHouseSpline
            error('loadCanadaData:housing', ...
                ['No monthly house-price series supplied.\n' ...
                 'FRED carries no monthly Canadian house-price index. Export the\n' ...
                 'Statistics Canada New Housing Price Index (table 18-10-0205, or\n' ...
                 'the Teranet-NBC HPI) to a two-column CSV and pass it as\n' ...
                 '   opts.houseFile = ''/path/to/nhpi_monthly.csv''\n' ...
                 'To run anyway on the spline-interpolated quarterly BIS proxy --\n' ...
                 'for code testing only, NOT for replication -- set\n' ...
                 '   opts.allowHouseSpline = true']);
        end
        % Any of the seven can be supplied from a file. opts.files is a
        % struct keyed by variable name; the older houseFile/cpiFile/rateFile
        % options map into it.
        userFile = extraFile;  userRow = extraRow;
        for j = 1:nV
            if isfield(opts.files, names{j}), userFile{j} = opts.files.(names{j}); end
            if isfield(opts.rows,  names{j}), userRow{j}  = opts.rows.(names{j});  end
        end
        legacy = {idx.house, opts.houseFile, opts.houseSeries; ...
                  idx.cpi,   opts.cpiFile,   opts.cpiSeries;   ...
                  idx.rate,  opts.rateFile,  opts.rateSeries};
        for a = 1:size(legacy,1)
            j = legacy{a,1};
            if isempty(userFile{j}) && ~isempty(legacy{a,2})
                userFile{j} = legacy{a,2};  userRow{j} = legacy{a,3};
            end
        end
        if isempty(userRow{idx.cpi})
            userRow{idx.cpi} = {'all-items','all items', ...
                                'ensemble','ensemble des biens et services'};
        end

        for j = 1:nV
            if ~isempty(userFile{j})
                uf = resolveData(opts.dataDir, userFile{j});   % look in data/ first
                fprintf('Reading %-16s (%s)\n', names{j}, uf);
                sj = readSeriesCSV(uf, names{j}, userRow{j});
                L(:,j) = placeSeries(sj, grid, false);
                ids{j} = uf;  quarterly(j) = false;
                fprintf('  %-6s %s .. %s, %d months read\n', names{j}, ...
                        datestr(sj.dates(1),'yyyy-mm'), datestr(sj.dates(end),'yyyy-mm'), ...
                        sum(~isnan(sj.value)));
            elseif strcmpi(source,'fred')
                fprintf('Fetching %-16s (%s)...\n', names{j}, ids{j});
                sj = fetchFRED(ids{j});
                L(:,j) = placeSeries(sj, grid, quarterly(j));
            else
                csvDir = opts.csvDir;
                if isempty(csvDir), csvDir = opts.dataDir; end   % fall back to data/
                if isempty(csvDir), csvDir = pwd; end
                f = fullfile(csvDir, [ids{j} '.csv']);
                if ~exist(f,'file')
                    error('loadCanadaData:local', ...
                        ['Missing %s.\nDownload it from ' ...
                         'https://fred.stlouisfed.org/graph/fredgraph.csv?id=%s ' ...
                         'and save it as %s.'], f, ids{j}, f);
                end
                fprintf('Reading %-16s (%s)\n', names{j}, f);
                sj = readFredCSV(f, ids{j});
                L(:,j) = placeSeries(sj, grid, quarterly(j));
            end
        end

        % Last observed month per series, so it is obvious which one caps the
        % panel and therefore which one to replace.
        fprintf('\n  last observation by series:\n');
        for j = 1:nV
            ok = find(~isnan(L(:,j)), 1, 'last');
            if isempty(ok)
                fprintf('    %-6s (none in range)\n', names{j});
            else
                fprintf('    %-6s %s\n', names{j}, datestr(grid(ok),'yyyy-mm'));
            end
        end
        fprintf('\n');

        if opts.allowHouseSpline && ~haveHouse
            fprintf(2, ['\n*** WARNING: housing is the spline-interpolated quarterly BIS proxy.\n' ...
                        '    Monthly growth from a splined quarterly index is ~95%% predictable\n' ...
                        '    from its own lags; Sigma_u, D, the IRFs, the variance decomposition\n' ...
                        '    and ALL scenarios are distorted. Results are not comparable to the\n' ...
                        '    paper. Supply opts.houseFile to fix. ***\n\n']);
        end

    otherwise
        error('loadCanadaData:source','Unknown source "%s".', source);
    end

    [data, span] = pack(grid, L, names, tcode, idx, true, opts.openEnded);

    if strcmpi(source,'fred') && ~isempty(opts.cacheFile)
        cachePath = resolveCache(opts.cacheDir, opts.cacheFile);
        cdir = fileparts(cachePath);
        if ~isempty(cdir) && ~exist(cdir,'dir'), mkdir(cdir); end
        save(cachePath,'data');
    end
end

% ======================================================================
function f = resolveData(dataDir, f)
% Resolve an INPUT data file. A bare filename ('nhpi.csv') is looked for in
% dataDir first; if found there it is used, otherwise it is left untouched so
% the old behaviour (relative to pwd) still works. Absolute paths pass through.
    if isempty(f) || isAbsolutePath(f) || isempty(dataDir), return; end
    cand = fullfile(dataDir, f);
    if exist(cand,'file'), f = cand; end
end

% ======================================================================
function f = resolveCache(cacheDir, f)
% Resolve a .mat CACHE file. A bare filename is placed under cacheDir (used
% for both load and save, so the file need not exist yet). A path that already
% carries a folder, or an absolute path, passes through unchanged.
    if isempty(f) || isAbsolutePath(f) || isempty(cacheDir), return; end
    if isempty(fileparts(f)), f = fullfile(cacheDir, f); end
end

% ======================================================================
function tf = isAbsolutePath(p)
    tf = false;
    if isempty(p), return; end
    if p(1) == '/' || p(1) == '\', tf = true; return; end          % unix / UNC
    if numel(p) >= 2 && isletter(p(1)) && p(2) == ':', tf = true; end  % windows drive
end

% ======================================================================
function [data, span] = pack(grid, L, names, tcode, idx, doClean, openEnded)
    if nargin < 7, openEnded = false; end
    span = seriesSpans(grid, L, names);
    if doClean
        [grid, L, span] = cleanPanel(grid, L, names, span, openEnded);
    end
    data = struct('dates',grid,'levels',L,'names',{names}, ...
                  'tcode',{tcode},'idx',idx);
end

function span = seriesSpans(grid, L, names)
    span = struct('name',{},'first',{},'last',{},'binding',{});
    for j = 1:size(L,2)
        ok = find(~isnan(L(:,j)));
        span(j).name    = names{j};
        span(j).binding = false;
        if isempty(ok)
            span(j).first = NaN;  span(j).last = NaN;
        else
            span(j).first = grid(ok(1));  span(j).last = grid(ok(end));
        end
    end
end

% ======================================================================
function [grid, L, span] = cleanPanel(grid, L, names, span, openEnded)
    if nargin < 5, openEnded = false; end
% Trim to the common non-NaN span, reporting exactly which series bind, then
% linearly fill sporadic interior gaps so OLS never sees NaNs.
    good = all(~isnan(L), 2);
    if ~any(good)
        reportSpans(span);
        error('loadCanadaData:empty', ...
            ['No month has all seven series present -- see the spans above. ' ...
             'At least one series does not overlap the others.']);
    end
    lo = find(good, 1, 'first');  hi = find(good, 1, 'last');

    starts = [span.first];  ends = [span.last];
    for j = 1:numel(span)
        span(j).binding = (starts(j) >= grid(lo)) || (ends(j) <= grid(hi));
    end

    if lo > 1 || hi < size(L,1)
        binding = names([span.last] == grid(hi));
        if openEnded
            % The trim is the point of an open-ended request, not a fault.
            fprintf('Panel capped at %s by: %s\n', ...
                    datestr(grid(hi),'yyyy-mm'), strjoin(binding, ', '));
        else
            fprintf(2, '\n*** Requested sample was TRIMMED: %s .. %s ***\n', ...
                    datestr(grid(lo),'yyyy-mm'), datestr(grid(hi),'yyyy-mm'));
            lateStart = find(starts > grid(1));
            earlyEnd  = find(ends   < grid(end));
            if ~isempty(lateStart)
                fprintf(2, '    starts late: %s\n', strjoin(names(lateStart), ', '));
            end
            if ~isempty(earlyEnd)
                fprintf(2, ['    ENDS EARLY : %s  <-- this moves the last observed month,\n' ...
                            '                 so every scenario would start from the wrong row.\n'], ...
                        strjoin(names(earlyEnd), ', '));
            end
            reportSpans(span);
        end
    end

    grid = grid(lo:hi);  L = L(lo:hi, :);

    for j = 1:size(L,2)
        col = L(:,j);  bad = isnan(col);
        if any(bad)
            fprintf(2,'  interior gaps in %s (%d months) -> linearly filled.\n', ...
                    names{j}, sum(bad));
            col(bad) = interp1(find(~bad), col(~bad), find(bad), 'linear', 'extrap');
            L(:,j) = col;
        end
    end
end

function reportSpans(span)
    fprintf('\n  %-8s %-10s %-10s\n', 'series', 'first', 'last');
    for j = 1:numel(span)
        if isnan(span(j).first)
            fprintf('  %-8s %-10s %-10s   (EMPTY)\n', span(j).name, '-', '-');
        else
            fprintf('  %-8s %-10s %-10s\n', span(j).name, ...
                    datestr(span(j).first,'yyyy-mm'), datestr(span(j).last,'yyyy-mm'));
        end
    end
    fprintf('\n');
end

% ======================================================================
function v = placeSeries(s, grid, isQuarterly)
    if isQuarterly
        v = quarterlyToMonthly(s.dates, s.value, grid);
    else
        v = alignMonthly(s.dates, s.value, grid);
    end
end

function v = alignMonthly(dates, value, grid)
    v = nan(numel(grid),1);
    [tf, loc] = ismember(grid, dates);
    v(tf) = value(loc(tf));
end

function v = quarterlyToMonthly(dates, value, grid)
    ok = ~isnan(value);
    xq = dates(ok);  yq = value(ok);
    inrng = grid >= min(xq) & grid <= max(xq);
    v = nan(numel(grid),1);
    v(inrng) = interp1(xq, yq, grid(inrng), 'spline');
end

% ======================================================================
function o = withDefaults(o, d)
    f = fieldnames(d);
    for i = 1:numel(f)
        if ~isfield(o, f{i})
            o.(f{i}) = d.(f{i});
        end
    end
end

% ======================================================================
function L = syntheticPanel(T)
% Plausible 7-variable monthly panel so the pipeline runs offline.
% NOT real data -- for code-testing and demonstration only.
    t = (1:T).';  randn('seed',42);
    oil   = max(40 + 30*sin(2*pi*t/120) + cumsum(randn(T,1))*1.5 + 40, 15);
    usip  = 100*exp(0.0015*t + 0.01*cumsum(randn(T,1)));
    cpi   = 100*exp(0.0018*t + 0.0008*cumsum(abs(randn(T,1))));
    house = 100*exp(0.0035*t + 0.004*cumsum(randn(T,1)));
    unemp = min(max(7 + 1.5*sin(2*pi*t/90) + cumsum(randn(T,1))/sqrt(T)*5,3),14);
    rate  = min(max(4 + 2*sin(2*pi*t/140) + cumsum(randn(T,1))/sqrt(T)*5,0.25),9);
    fx    = 1.3 + 0.15*sin(2*pi*t/100) + cumsum(randn(T,1))/sqrt(T)*3*0.02;
    L = [oil usip cpi house unemp rate fx];
end
