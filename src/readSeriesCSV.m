function S = readSeriesCSV(path, label, rowLabel)
% READSERIESCSV  Read a monthly series from CSV into {dates(datenum), value}.
%
%   S = readSeriesCSV('nhpi.csv')
%   S = readSeriesCSV('1810020501-eng.csv', 'house')
%   S = readSeriesCSV('1810020501-eng.csv', 'house', 'House only')
%
% Two layouts are auto-detected.
%
% LONG -- one row per month, two or more columns:
%       observation_date,VALUE
%       1992-01-01,55.8
%   Dates may be 'yyyy-mm-dd', 'yyyy-mm', 'yyyy/mm/dd' or 'mm/dd/yyyy';
%   missing values may be '.', '..', 'NA', 'NaN', 'F', 'x' or empty.
%
%   With MORE than two columns -- a Bank of Canada Valet export, say, whose
%   observation block carries every sub-index side by side -- name the column
%   you want in the third argument:
%
%       readSeriesCSV('bcpi.csv', 'oil', 'M.ENER')
%
%   Omitting it errors with the list of columns rather than defaulting to the
%   second one, because silently returning the wrong sub-index is worse than
%   failing. Any preamble above the observation block is skipped.
%
% WIDE -- Statistics Canada's "Download as displayed" export, in which the
% months run ACROSS a header row and each series is one row beneath it:
%       "New housing price indexes","January 1981","February 1981",...
%       "Total (house and land)","38.2","38.7",...
%       "House only","36.1","36.5",...
%   The leading metadata block, the trailing symbol legend, table
%   corrections and footnotes are all skipped. Data-quality flags appended
%   to a value ('40.6E') are stripped. English and French month names are
%   both recognised, so the -eng and -fra downloads work interchangeably.
%
% INPUTS
%   path      CSV file
%   label     name recorded in S.id (default: the path)
%   rowLabel  which series to take. In the WIDE layout it names a ROW
%             (default 'Total (house and land)'); in a multi-column LONG
%             layout it names a COLUMN. Either way, a label that is not found
%             produces an error listing what the file actually contains.
%
% OUTPUT
%   S.id, S.dates (datenum, normalised to the 1st of the month), S.value

    if nargin < 2 || isempty(label),    label = path; end
    if nargin < 3,                      rowLabel = ''; end

    lines = readLines(path);
    cells = cell(numel(lines),1);
    for i = 1:numel(lines)
        cells{i} = splitCSVLine(lines{i});
    end

    hdrRow = findWideHeader(cells);
    if isempty(hdrRow)
        [dates, value] = parseLong(cells, path, rowLabel);
    else
        [dates, value] = parseWide(cells, hdrRow, rowLabel, path);
    end

    [dates, ord] = sort(dates);  value = value(ord);
    S = struct('id', label, 'dates', dates, 'value', value);
end

% ======================================================================
%                              LONG LAYOUT
% ======================================================================
function [dates, value] = parseLong(cells, path, colLabel)
    n = numel(cells);

    % First data row, and the widest data row -- a preamble may hold narrow
    % rows that happen to start with something date-like.
    first = 0;  nCol = 0;
    for i = 1:n
        f = cells{i};
        if numel(f) >= 2 && ~isnan(parseDate(f{1}))
            if first == 0, first = i; end
            nCol = max(nCol, numel(f));
        end
    end
    if first == 0
        error('readSeriesCSV:parse', ...
            ['No parseable rows in %s, and it does not look like a ' ...
             'Statistics Canada wide export either. Expected a monthly date ' ...
             'and at least one numeric column.'], path);
    end

    % Column names come from the last non-data row above the block.
    colNames = {};
    if first > 1
        h = cells{first-1};
        if numel(h) >= 2 && isnan(parseDate(h{1}))
            colNames = cellfun(@(x) strtrim(strrep(x,'"','')), h, ...
                               'UniformOutput', false);
        end
    end

    % Which column?
    col = 2;
    if nCol > 2
        if isempty(colLabel)
            if isempty(colNames)
                error('readSeriesCSV:longAmbiguous', ...
                    ['%s has %d value columns and no header to name them. ' ...
                     'Pass a column index as the third argument.'], path, nCol-1);
            end
            error('readSeriesCSV:longAmbiguous', ...
                ['%s has %d value columns; name the one you want as the third ' ...
                 'argument.\nColumns available: %s'], path, nCol-1, ...
                 strjoin(colNames(2:min(nCol,numel(colNames))), ' | '));
        end
        if isnumeric(colLabel)
            col = colLabel;
        else
            if iscell(colLabel), want = colLabel; else, want = {colLabel}; end
            col = 0;
            for c = 2:numel(colNames)
                if anyMatch(colNames{c}, want, true), col = c; break; end
            end
            if col == 0
                for c = 2:numel(colNames)
                    if anyMatch(colNames{c}, want, false), col = c; break; end
                end
            end
            if col == 0
                error('readSeriesCSV:longColNotFound', ...
                    'Column "%s" not found in %s.\nColumns available: %s', ...
                    strjoin(want,'" / "'), path, ...
                    strjoin(colNames(2:min(nCol,numel(colNames))), ' | '));
            end
        end
    end

    dates = zeros(n,1);  value = nan(n,1);  keep = false(n,1);
    for i = first:n
        f = cells{i};
        if numel(f) < col, continue; end
        dn = parseDate(f{1});
        if isnan(dn), continue; end
        dates(i) = dn;  value(i) = parseValue(f{col});  keep(i) = true;
    end
    dates = dates(keep);  value = value(keep);
end

% ======================================================================
%                              WIDE LAYOUT
% ======================================================================
function hdrRow = findWideHeader(cells)
% The date header is the first row with many fields that parse as month-year.
    hdrRow = [];
    for i = 1:numel(cells)
        f = cells{i};
        if numel(f) < 13, continue; end
        nd = 0;
        for j = 2:numel(f)
            if ~isnan(parseMonthYear(f{j})), nd = nd + 1; end
        end
        if nd >= 12 && nd >= 0.8*(numel(f)-1)
            hdrRow = i;  return;
        end
    end
end

function [dates, value] = parseWide(cells, hdrRow, rowLabel, path)
    hdr = cells{hdrRow};
    col = [];  dts = [];
    for j = 2:numel(hdr)
        dn = parseMonthYear(hdr{j});
        if ~isnan(dn), col(end+1) = j; dts(end+1) = dn; end %#ok<AGROW>
    end

    if isempty(rowLabel)
        wanted = {'total (house and land)','total (maison et terrain)'};
    elseif iscell(rowLabel)
        wanted = cell(1,numel(rowLabel));
        for w = 1:numel(rowLabel), wanted{w} = lower(strtrim(rowLabel{w})); end
    else
        wanted = {lower(strtrim(rowLabel))};
    end

    % Rows that carry data: a non-empty first field and mostly numeric cells
    % in the date columns. This excludes the unit row, the symbol legend,
    % the table-corrections block and the footnotes.
    labels = {};  vals = {};
    for i = hdrRow+1:numel(cells)
        f = cells{i};
        if numel(f) < max(col), continue; end
        lab = strtrim(f{1});
        if isempty(lab), continue; end
        v = nan(numel(col),1);
        for k = 1:numel(col), v(k) = parseValue(f{col(k)}); end
        if mean(~isnan(v)) < 0.5, continue; end
        labels{end+1} = lab;  vals{end+1} = v; %#ok<AGROW>
    end

    if isempty(labels)
        error('readSeriesCSV:wideNoData', ...
            'Found a month header in %s but no data rows beneath it.', path);
    end

    % Exact match on the normalised label wins. Only if nothing matches
    % exactly do we fall back to a prefix match, and an ambiguous prefix is
    % an error rather than a silent choice -- asking for 'All-items' must not
    % quietly return 'All-items excluding food and energy'.
    hit = 0;
    for i = 1:numel(labels)
        if anyMatch(labels{i}, wanted, true), hit = i; break; end
    end
    if hit == 0
        cand = [];
        for i = 1:numel(labels)
            if anyMatch(labels{i}, wanted, false), cand(end+1) = i; end %#ok<AGROW>
        end
        if numel(cand) == 1
            hit = cand;
        elseif numel(cand) > 1
            error('readSeriesCSV:wideRowAmbiguous', ...
                ['Row "%s" matches %d rows in %s: %s\n' ...
                 'Give the full label instead.'], strjoin(wanted, '" / "'), ...
                 numel(cand), path, strjoin(labels(cand), ' | '));
        end
    end

    if hit == 0
        error('readSeriesCSV:wideRowNotFound', ...
            ['Row "%s" not found in %s.\nRows available: %s\n' ...
             'Pass the one you want as the third argument to readSeriesCSV ' ...
             '(opts.houseSeries / opts.cpiSeries).'], ...
             strjoin(wanted, '" / "'), path, strjoin(labels, ' | '));
    end

    dates = dts(:);  value = vals{hit}(:);
end

function tf = anyMatch(lab, wanted, exact)
% Compare a row label against the candidates, ignoring case, surrounding
% whitespace and any trailing footnote reference numbers. Statistics Canada
% appends those in the "as displayed" export, so the CPI all-items row
% arrives as 'All-items 8' and its header as 'Products and product groups 7'.
    tf = false;
    L = normLabel(lab);
    for k = 1:numel(wanted)
        W = normLabel(wanted{k});
        if isempty(W), continue; end
        if exact
            if strcmp(L, W), tf = true; return; end
        else
            if numel(L) >= numel(W) && strcmp(L(1:numel(W)), W), tf = true; return; end
        end
    end
end

function t = normLabel(s)
    t = lower(strtrim(s));
    stripped = regexprep(t, '(\s+\d+)+$', '');   % drop trailing footnote marks
    stripped = strtrim(stripped);
    if ~isempty(stripped) && any(isletter(stripped))
        t = stripped;                            % keep labels that are only digits
    end
end

% ======================================================================
%                              PRIMITIVES
% ======================================================================
function lines = readLines(path)
    fid = fopen(path, 'r');
    if fid < 0, error('readSeriesCSV:open', 'Cannot open %s', path); end
    raw = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', ''); fclose(fid);
    lines = raw{1};
    if ~isempty(lines)
        b = double(lines{1});                       % strip a UTF-8 BOM
        if numel(b) >= 3 && all(b(1:3) == [239 187 191])
            lines{1} = char(b(4:end));
        end
    end
    for i = 1:numel(lines)
        s = lines{i};
        s(double(s) == 160) = ' ';                  % non-breaking -> plain
        while ~isempty(s) && (s(end) == char(13)), s(end) = []; end   % CRLF
        lines{i} = s;
    end
    lines = lines(~cellfun(@isempty, lines));
end

function f = splitCSVLine(line)
% Split on commas that are not inside double quotes; unescape "" -> ".
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

function v = parseValue(s)
% Numeric value with any trailing data-quality flag removed ('40.6E' -> 40.6).
    v = NaN;
    s = strtrim(s);
    if isempty(s), return; end
    s = s(s ~= ',');                    % StatCan writes levels as "16,852.4"
    while ~isempty(s) && isletter(s(end)), s(end) = []; end
    s = strtrim(s);
    if isempty(s) || all(s == '.'), return; end
    d = str2double(s);
    if ~isnan(d), v = d; end
end

function dn = parseDate(ds)
% First-of-month datenum for a LONG-layout date field, or NaN.
    dn = NaN;
    ds = strtrim(ds);
    if isempty(ds), return; end
    fmts = {'yyyy-mm-dd','yyyy-mm','yyyy/mm/dd','mm/dd/yyyy'};
    for k = 1:numel(fmts)
        try
            dv = datevec(ds, fmts{k});
            if dv(1) > 1800 && dv(1) < 2200 && dv(2) >= 1 && dv(2) <= 12
                dn = datenum(dv(1), dv(2), 1);  return;
            end
        catch
            % try the next format
        end
    end
    dn = parseMonthYear(ds);            % also accept 'January 1981'
end

function dn = parseMonthYear(s)
% First-of-month datenum for 'January 1981' / 'janvier 1981', or NaN.
% Parsed by hand rather than through datevec's 'mmmm', whose locale
% behaviour differs between MATLAB and Octave.
    persistent MON
    if isempty(MON)
        MON = { ...
          {'january','jan','janvier','janv'}, ...
          {'february','feb','fevrier','février','fevr','févr'}, ...
          {'march','mar','mars'}, ...
          {'april','apr','avril','avr'}, ...
          {'may','mai'}, ...
          {'june','jun','juin'}, ...
          {'july','jul','juillet','juil'}, ...
          {'august','aug','aout','août'}, ...
          {'september','sep','sept','septembre'}, ...
          {'october','oct','octobre'}, ...
          {'november','nov','novembre'}, ...
          {'december','dec','decembre','décembre'} };
    end
    dn = NaN;
    s = lower(strtrim(s));
    if isempty(s), return; end
    tok = regexp(s, '^([a-zA-Zàâçéèêëîïôûùüÿñæœ\.]+)\s+(\d{4})$', 'tokens', 'once');
    if isempty(tok), return; end
    name = strrep(tok{1}, '.', '');
    yr   = str2double(tok{2});
    for m = 1:12
        if any(strcmp(name, MON{m}))
            dn = datenum(yr, m, 1);  return;
        end
    end
end
