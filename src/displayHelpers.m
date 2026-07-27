function h = displayHelpers()
% DISPLAYHELPERS  Returns a struct of function handles for figure display
% conversions. Keeps main_risk_scenarios readable.
%
%   d = displayHelpers();  d.toLevels(...);  d.toYoY(...);  d.annualize(...)

    h.toLevels  = @toLevels;
    h.toYoY     = @toYoY;
    h.annualize = @annualize;
end

function lev = toLevels(lastLevel, dlogPath)
% Reconstruct an index level from a starting level and monthly log-diffs.
    lev = lastLevel * exp(cumsum(dlogPath(:)));
end

function y = annualize(dlogPath)
% Annualised monthly growth in percent: (e^{1200*mean}?) -- simple *1200.
    y = 1200*dlogPath(:);
end

function yoy = toYoY(levelSeries)
% 12-month percent change of a level series (NaN for first 12 points).
    n = numel(levelSeries); yoy = nan(n,1);
    yoy(13:end) = 100*(levelSeries(13:end)./levelSeries(1:end-12) - 1);
end
