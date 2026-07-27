function [dax, dval] = dispSeries(data, v, mode, fc, H)
% DISPSERIES  Assemble a history+forecast display series for variable v.
%
%   [dax, dval] = dispSeries(data, v, mode, fc, H)
%
% INPUTS
%   data  raw data struct (dates datenum, levels, tcode)
%   v     variable column index
%   mode  'level' | 'yoy' | 'agrowth'
%   fc    H x 1 forecast in transformed units:
%           - for dlog variables: forecast monthly log-differences
%           - for level variables: forecast levels
%   H     horizon
%
% OUTPUT
%   dax   datenum axis (history then forecast)
%   dval  display values on that axis

    rawd = data.dates;  L = data.levels;  tc = data.tcode{v};
    fdates = zeros(H,1); dv = datevec(rawd(end));
    for h = 1:H, fdates(h) = datenum(dv(1), dv(2)+h, 1); end

    switch tc
    case 'dlog'
        flev = L(end,v) * exp(cumsum(fc(:)));         % reconstructed levels
        switch mode
        case 'level'
            dax = [rawd; fdates];        dval = [L(:,v); flev];
        case 'yoy'
            comb = [L(:,v); flev];
            dax  = [rawd; fdates];       dval = yoy12(comb);
        case 'agrowth'
            gh   = 1200*diff(log(L(:,v)));            % history, aligned rawd(2:end)
            dax  = [rawd(2:end); fdates]; dval = [gh; 1200*fc(:)];
        otherwise, error('dispSeries:mode','bad mode');
        end
    case 'level'
        dax = [rawd; fdates];            dval = [L(:,v); fc(:)];
    end
end

function y = yoy12(levelSeries)
    n = numel(levelSeries); y = nan(n,1);
    y(13:end) = 100*(levelSeries(13:end)./levelSeries(1:end-12) - 1);
end
