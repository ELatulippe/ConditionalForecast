function [X, labels] = covidDummies(dates, spec)
% COVIDDUMMIES  Month dummies for the pandemic episode, as VAR regressors.
%
%   X = covidDummies(dates)                  % no dummies (spec empty)
%   X = covidDummies(dates, 'default')       % March-August 2020
%   X = covidDummies(dates, [2020 3; 2020 8])
%
% WHY THIS EXISTS. In April-May 2020 the oil price, U.S. industrial
% production and the Canadian unemployment rate all moved by many times their
% historical standard deviations, in the same months. Those observations
% dominate Sigma_u. Because the recursive identification orders oil first and
% U.S. activity second, the common pandemic movement is attributed to those
% two shocks, which inflates their share of the unemployment forecast-error
% variance and pushes unemployment's own share to the floor. The same episode
% flattens the estimated oil-to-CPI pass-through, since a 57% monthly oil
% innovation is not matched proportionally in consumer prices.
%
% One dummy per month over the specified window absorbs those observations,
% so Sigma_u and the impulse responses are estimated off the rest of the
% sample. The dummies are zero outside the window and therefore zero over the
% forecast horizon, which leaves the baseline forecast recursion untouched.
%
% THIS IS A DIAGNOSTIC, NOT A REPLICATION. The paper does not mention
% pandemic dummies. Use this to find out whether the 2020 episode is what
% separates your estimates from Table 1; do not present dummied results as a
% reproduction of the paper.
%
% INPUTS
%   dates  T x 1 datenums, one per row of the VAR's data matrix
%   spec   [] or '' for none; 'default' for March-August 2020; or a
%          2 x 2 matrix [startYear startMonth; endYear endMonth]
%
% OUTPUTS
%   X       T x q dummy matrix (empty when spec is empty)
%   labels  1 x q cellstr, 'yyyy-mm' for each dummy

    X = [];  labels = {};
    if nargin < 2 || isempty(spec), return; end

    if ischar(spec)
        if strcmpi(spec, 'none'), return; end
        if ~strcmpi(spec, 'default')
            error('covidDummies:spec', ...
                'Unknown spec "%s". Use ''default'', ''none'', or [y1 m1; y2 m2].', spec);
        end
        spec = [2020 3; 2020 8];
    end
    if ~isequal(size(spec), [2 2])
        error('covidDummies:spec', ...
            'Numeric spec must be 2x2: [startYear startMonth; endYear endMonth].');
    end

    lo = datenum(spec(1,1), spec(1,2), 1);
    hi = datenum(spec(2,1), spec(2,2), 1);
    if hi < lo
        error('covidDummies:spec', 'End of the window precedes its start.');
    end

    sel = find(dates >= lo & dates <= hi);
    if isempty(sel)
        warning('covidDummies:outOfSample', ...
            ['The window %s..%s does not overlap the estimation sample; ' ...
             'no dummies added.'], datestr(lo,'yyyy-mm'), datestr(hi,'yyyy-mm'));
        return;
    end

    T = numel(dates);
    X = zeros(T, numel(sel));
    labels = cell(1, numel(sel));
    for k = 1:numel(sel)
        X(sel(k), k) = 1;
        labels{k} = datestr(dates(sel(k)), 'yyyy-mm');
    end
end
