function [L, ix, originMonth, R0] = scenarioOrigin(ob)
%SCENARIOORIGIN  Pass 1 -- read the forecast origin with no bands, no figures.
%
%   [L, ix, originMonth, R0] = scenarioOrigin(ob)
%
% Every scenario reads its starting level off the last row of the panel
% (oil0, u0, r0, ...), so before building any conditioning path we estimate the
% VAR once, cheaply, to recover those terminal levels and the variable index.
% This is the "pass 1" step from Running_code.txt, factored out so that the
% single / joint / temporary / soft / counterfactual / FEVD scripts can each
% call it instead of repeating the block.
%
% INPUT
%   ob   options struct, typically from scenario_config()
%
% OUTPUT
%   L            1 x K row of terminal LEVELS (index with ix.<name>)
%   ix           struct mapping variable name -> column, e.g. ix.oil
%   originMonth  datenum of the forecast origin (last month of the panel)
%   R0           the full return struct from main_risk_scenarios (has .model,
%                .tr, .data) -- handy for counterfactuals and the FEVD

    ob.bands   = false;   % pass 1: no uncertainty bands ...
    ob.figures = false;   % ... and no figures, just the estimate

    R0 = main_risk_scenarios('fred', ob);

    L           = R0.data.levels(end,:);
    ix          = R0.data.idx;
    originMonth = R0.data.dates(end);

    fprintf('origin %s: oil %.1f, unemp %.2f%%, rate %.2f%%\n', ...
            datestr(originMonth,'yyyy-mm'), L(ix.oil), L(ix.unemp), L(ix.rate));
end
