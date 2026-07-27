function S = fetchFRED(seriesID, opts)
% FETCHFRED  Download one FRED series as a CSV (no API key required).
%
%   S = fetchFRED('EXCAUS')
%   S = fetchFRED('EXCAUS', struct('timeout',60,'retries',5))
%
% Uses the public fredgraph CSV endpoint:
%   https://fred.stlouisfed.org/graph/fredgraph.csv?id=<ID>
%
% OPTS (all optional):
%   .timeout  seconds per attempt   (default 60)
%   .retries  number of attempts    (default 4)
%
% OUTPUT struct S with .id, .dates (datenum, 1st of month), .value (NaN='.').
%
% If your network blocks or throttles FRED, download the CSVs once in a
% browser and use loadCanadaData(...,'local', csvDir) instead (see README).

    if nargin < 2, opts = struct(); end
    if ~isfield(opts,'timeout'), opts.timeout = 60; end
    if ~isfield(opts,'retries'), opts.retries = 4;  end

    url = ['https://fred.stlouisfed.org/graph/fredgraph.csv?id=' seriesID];
    tmp = [tempname '.csv'];

    % Longer timeout + a browser-like user agent avoid the default 5s cut-off
    % and the occasional bot filter.
    wo = weboptions('Timeout', opts.timeout, 'ContentType', 'text', ...
                    'UserAgent', 'Mozilla/5.0 (MATLAB fetchFRED)');

    lastErr = '';
    for attempt = 1:opts.retries
        try
            websave(tmp, url, wo);
            S = readFredCSV(tmp, seriesID);
            if exist(tmp,'file'), delete(tmp); end
            return;                                   % success
        catch ME
            lastErr = ME.message;
            if exist(tmp,'file'), delete(tmp); end
            pause(2*attempt);                         % back off, then retry
        end
    end

    error('fetchFRED:download', ...
        ['Could not download %s after %d attempts (%s).\n' ...
         'Fixes: (1) rerun -- FRED is often slow on the first hit; ' ...
         '(2) raise the timeout, e.g. fetchFRED(''%s'',struct(''timeout'',120)); ' ...
         '(3) if FRED is blocked on your network, download the CSVs in a ' ...
         'browser and use loadCanadaData(...,''local'',csvDir).'], ...
         seriesID, opts.retries, lastErr, seriesID);
end
