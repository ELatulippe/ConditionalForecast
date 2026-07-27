function outFile = fetchStatCan(tableId, filters, outFile, opts)
% FETCHSTATCAN  Download a Statistics Canada table and extract one series.
%
%   fetchStatCan('10100144')                                  % list what's in it
%   fetchStatCan('10100144', {'Rates','Bank rate'}, 'rate.csv')
%   fetchStatCan('14100287', {'Labour force characteristics','Unemployment rate'; ...
%                             'Data type','Seasonally adjusted'}, 'unemp.csv')
%
% Writes a two-column date,value CSV that plugs straight into
%   opts.files = struct('rate','rate.csv')
%
% StatCan publishes every table as a zipped full-table CSV at
%   https://www150.statcan.gc.ca/n1/tbl/csv/<8-digit id>-eng.zip
% in LONG format: one row per period per series, with the series identified by
% label columns (GEO, the table's dimension columns, and so on). This function
% downloads it, filters the rows, and keeps REF_DATE and VALUE.
%
% Called with fewer than three arguments, it downloads and then PRINTS the
% distinct values of every label column instead of writing anything, so you
% can see what to filter on. Start there -- dimension names differ by table.
%
% INPUTS
%   tableId  '10100144', '10-10-0144-01' or '1010014401' (all accepted)
%   filters  Nx2 cell array of {columnName, value} pairs, matched
%            case-insensitively after trimming and after dropping any
%            trailing footnote digits StatCan appends to labels
%   outFile  path of the two-column CSV to write
%   opts     .localZip  use this .zip instead of downloading
%            .localCsv  use this .csv instead of downloading (skips the unzip)
%            .geo       shorthand for a {'GEO', value} filter (default 'Canada')
%            .timeout   seconds for the download (default 120)
%
% NOTE ON UNITS. Rate tables are usually already in percent, matching what the
% VAR expects for the policy rate. Check the UOM column in the listing before
% using anything unfamiliar.

    if nargin < 2, filters = {}; end
    if nargin < 3, outFile = ''; end
    if nargin < 4, opts = struct(); end
    if ~isfield(opts,'localZip'), opts.localZip = ''; end
    if ~isfield(opts,'localCsv'), opts.localCsv = ''; end
    if ~isfield(opts,'geo'),      opts.geo      = 'Canada'; end
    if ~isfield(opts,'timeout'),  opts.timeout  = 120; end

    id8 = normaliseId(tableId);

    % ---- obtain the CSV -------------------------------------------------
    cleanupDir = '';
    if ~isempty(opts.localCsv)
        csvPath = opts.localCsv;
    else
        if ~isempty(opts.localZip)
            zipPath = opts.localZip;
        else
            url = sprintf('https://www150.statcan.gc.ca/n1/tbl/csv/%s-eng.zip', id8);
            zipPath = [tempname '.zip'];
            fprintf('Downloading %s ...\n', url);
            try
                wo = weboptions('Timeout', opts.timeout, ...
                                'UserAgent','Mozilla/5.0 (MATLAB fetchStatCan)');
                websave(zipPath, url, wo);
            catch ME
                error('fetchStatCan:download', ...
                    ['Could not download table %s (%s).\n' ...
                     'Download it in a browser from\n  %s\n' ...
                     'then pass opts.localZip or opts.localCsv.'], ...
                     id8, ME.message, url);
            end
        end
        cleanupDir = [tempname '_sc'];
        mkdir(cleanupDir);
        unzip(zipPath, cleanupDir);
        d = dir(fullfile(cleanupDir, '*.csv'));
        keep = ~cellfun(@(n) ~isempty(strfind(lower(n),'metadata')), {d.name});
        d = d(keep);
        if isempty(d)
            error('fetchStatCan:nocsv', 'No data CSV inside the archive.');
        end
        [~, big] = max([d.bytes]);
        csvPath = fullfile(cleanupDir, d(big).name);
    end

    % ---- parse ----------------------------------------------------------
    [hdr, rows] = readLongCsv(csvPath);
    cDate  = findCol(hdr, 'REF_DATE');
    cValue = findCol(hdr, 'VALUE');
    if isempty(cDate) || isempty(cValue)
        error('fetchStatCan:cols', ...
            'Expected REF_DATE and VALUE columns; found: %s', strjoin(hdr, ', '));
    end

    % Label columns = everything that is not a date, a value, or bookkeeping.
    skip = {'REF_DATE','VALUE','DGUID','UOM_ID','SCALAR_ID','VECTOR', ...
            'COORDINATE','STATUS','SYMBOL','TERMINATED','DECIMALS','SCALAR_FACTOR'};
    isLabel = true(1, numel(hdr));
    for c = 1:numel(hdr)
        if any(strcmpi(hdr{c}, skip)), isLabel(c) = false; end
    end

    % ---- listing mode ---------------------------------------------------
    if isempty(outFile)
        fprintf('\nTable %s -- distinct values per label column:\n', id8);
        for c = find(isLabel)
            vals = uniqueStrings(rows(:,c));
            fprintf('\n  %s  (%d values)\n', hdr{c}, numel(vals));
            for k = 1:min(numel(vals), 25)
                fprintf('    %s\n', vals{k});
            end
            if numel(vals) > 25, fprintf('    ... and %d more\n', numel(vals)-25); end
        end
        fprintf('\nPass the ones you want as filters, e.g.\n');
        fprintf('  fetchStatCan(''%s'', {''%s'',''<value>''}, ''out.csv'')\n\n', ...
                id8, hdr{find(isLabel,1,'last')});
        if ~isempty(cleanupDir), rmdir(cleanupDir,'s'); end
        return;
    end

    % ---- filter ---------------------------------------------------------
    if ~isempty(opts.geo)
        cGeo = findCol(hdr, 'GEO');
        if ~isempty(cGeo)
            filters = [{'GEO', opts.geo}; filters];
        end
    end

    keep = true(size(rows,1), 1);
    for f = 1:size(filters,1)
        c = findCol(hdr, filters{f,1});
        if isempty(c)
            error('fetchStatCan:filtercol', ...
                'No column "%s". Columns: %s', filters{f,1}, strjoin(hdr, ', '));
        end
        want = normLab(filters{f,2});
        hit  = false(size(rows,1),1);
        for r = 1:size(rows,1)
            hit(r) = strcmp(normLab(rows{r,c}), want);
        end
        if ~any(hit)
            vals = uniqueStrings(rows(keep,c));
            error('fetchStatCan:nomatch', ...
                'No row has %s = "%s". Values present: %s', ...
                filters{f,1}, filters{f,2}, strjoin(vals(1:min(end,20)), ' | '));
        end
        keep = keep & hit;
    end
    if ~any(keep)
        error('fetchStatCan:empty', 'Filters matched no rows jointly.');
    end

    dates = rows(keep, cDate);
    vals  = rows(keep, cValue);
    if numel(dates) ~= numel(unique(dates))
        left = setdiff(hdr(isLabel), [{'GEO'}, filters(:,1)']);
        error('fetchStatCan:ambiguous', ...
            ['The filters leave %d rows for %d distinct dates, so more than ' ...
             'one series is still selected. Add a filter on: %s'], ...
             numel(dates), numel(unique(dates)), strjoin(left, ', '));
    end

    % ---- write ----------------------------------------------------------
    fid = fopen(outFile, 'w');
    if fid < 0, error('fetchStatCan:write', 'Cannot write %s', outFile); end
    fprintf(fid, 'date,value\n');
    for r = 1:numel(dates)
        fprintf(fid, '%s,%s\n', strtrim(dates{r}), strtrim(vals{r}));
    end
    fclose(fid);
    fprintf('Wrote %s: %d observations, %s .. %s\n', ...
            outFile, numel(dates), strtrim(dates{1}), strtrim(dates{end}));

    if ~isempty(cleanupDir), rmdir(cleanupDir,'s'); end
end

% ======================================================================
function id8 = normaliseId(t)
    d = t(t >= '0' & t <= '9');
    if numel(d) < 8
        error('fetchStatCan:id', 'Table id "%s" does not contain 8 digits.', t);
    end
    id8 = d(1:8);
end

function c = findCol(hdr, name)
    c = find(strcmpi(strtrim(hdr), name), 1);
end

function t = normLab(s)
    t = lower(strtrim(s));
    t = t(t ~= '"');
    stripped = strtrim(regexprep(t, '(\s+\d+)+$', ''));
    if ~isempty(stripped) && any(isletter(stripped)), t = stripped; end
end

function u = uniqueStrings(col)
    n = cellfun(@normLab, col, 'UniformOutput', false);
    [~, ia] = unique(n);
    u = col(sort(ia));
    u = cellfun(@(x) strtrim(x), u, 'UniformOutput', false);
end

% ======================================================================
function [hdr, rows] = readLongCsv(path)
    fid = fopen(path, 'r');
    if fid < 0, error('fetchStatCan:open', 'Cannot open %s', path); end
    raw = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', ''); fclose(fid);
    lines = raw{1};
    lines = lines(~cellfun(@isempty, lines));
    if isempty(lines), error('fetchStatCan:empty', '%s is empty.', path); end

    b = double(lines{1});
    if numel(b) >= 3 && all(b(1:3) == [239 187 191]), lines{1} = char(b(4:end)); end

    hdr  = splitCSVLine(lines{1});
    hdr  = cellfun(@(x) strtrim(strrep(x,'"','')), hdr, 'UniformOutput', false);
    nCol = numel(hdr);
    rows = cell(numel(lines)-1, nCol);
    keep = false(numel(lines)-1, 1);
    for i = 2:numel(lines)
        f = splitCSVLine(lines{i});
        if numel(f) < nCol, continue; end
        rows(i-1, :) = f(1:nCol);
        keep(i-1) = true;
    end
    rows = rows(keep, :);
end

function f = splitCSVLine(line)
    f = {};  cur = '';  inq = false;  i = 1;  n = numel(line);
    while i <= n
        c = line(i);
        if c == '"'
            if inq && i < n && line(i+1) == '"'
                cur(end+1) = '"'; i = i + 2; continue; %#ok<AGROW>
            end
            inq = ~inq;  i = i + 1;  continue;
        end
        if c == ',' && ~inq
            f{end+1} = cur; cur = ''; i = i + 1; continue; %#ok<AGROW>
        end
        cur(end+1) = c;  i = i + 1; %#ok<AGROW>
    end
    f{end+1} = cur;
end
