function tr = transformData(data)
% TRANSFORMDATA  Map raw levels to the stationary variables the VAR uses.
%
%   tr = transformData(data)
%
% First difference of logs for oil, usip, cpi, house, fx; levels for unemp
% and rate (matching the paper). Monthly log-differences are stored raw
% (NOT annualised); annualisation is a display choice handled elsewhere.
%
% OUTPUT struct 'tr':
%   .dates   (T-1) x 1 datetime aligned to the transformed rows
%   .Y       (T-1) x 7 stationary data (VAR input)
%   .names,.tcode,.idx  copied from data
%   .lastLevels 1 x 7   final observed raw levels (for level reconstruction)
%   .levels  T x 7 raw levels (kept for scenario construction)
%   .rawdates T x 1

    L = data.levels;  tcode = data.tcode;  [T,K] = size(L);
    Y = zeros(T-1, K);
    for j = 1:K
        switch tcode{j}
            case 'dlog'
                Y(:,j) = diff(log(L(:,j)));
            case 'level'
                Y(:,j) = L(2:end, j);
            otherwise
                error('transformData:tcode','Unknown tcode %s', tcode{j});
        end
    end
    tr = struct('dates', data.dates(2:end), 'Y', Y, ...
                'names',{data.names}, 'tcode',{data.tcode}, 'idx',data.idx, ...
                'lastLevels', L(end,:), 'levels', L, 'rawdates', data.dates);
end
