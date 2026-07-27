function S = readFredCSV(path, seriesID)
% READFREDCSV  Parse a FRED CSV file into {dates(datenum), value}.
% Handles both the legacy "DATE,<ID>" and newer "observation_date,<ID>"
% headers, and missing values encoded as '.'.
%
%   S = readFredCSV('MCOILWTICO.csv', 'MCOILWTICO')

    fid = fopen(path,'r');
    if fid < 0, error('readFredCSV:open','Cannot open %s', path); end
    raw = textscan(fid,'%s','Delimiter','\n'); fclose(fid);
    lines = raw{1};
    lines = lines(~cellfun(@isempty,lines));
    n = numel(lines) - 1;
    dates = zeros(n,1);  value = nan(n,1);
    for i = 2:numel(lines)                       % skip header row
        parts = strsplit(lines{i}, ',');
        if numel(parts) < 2, continue; end       % skip malformed rows
        ds = strtrim(parts{1});
        vs = strtrim(parts{2});
        dv = datevec(ds, 'yyyy-mm-dd');
        dates(i-1) = datenum(dv(1), dv(2), 1);   % normalise to 1st of month
        if ~strcmp(vs,'.') && ~isempty(vs)
            value(i-1) = str2double(vs);
        end
    end
    S = struct('id', seriesID, 'dates', dates, 'value', value);
end
