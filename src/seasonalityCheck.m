function S = seasonalityCheck(tr, verbose)
% SEASONALITYCHECK  Test each growth-rate series for deterministic monthly
% seasonality, by regressing it on a constant and eleven month dummies.
%
%   S = seasonalityCheck(tr)
%   S = seasonalityCheck(tr, false)   % silent
%
% Why this matters here. The paper's database is seasonally adjusted, but the
% FRED default for Canadian CPI (CPALCY01CAM661N) is not. A deterministic
% seasonal in a monthly log-difference is, from the VAR's point of view, a
% large component of that variable's innovation that is orthogonal to every
% other variable. It therefore inflates the variable's OWN share of its
% forecast-error variance and deflates everything else's -- which is exactly
% the signature to look for when a row of Table 1 is too self-explained.
%
% A twelve-month growth rate nets deterministic seasonality out, so the
% figures can look entirely reasonable while the estimated monthly dynamics
% are wrong. This check runs on the estimation-frequency data instead.
%
% INPUT
%   tr       struct from transformData
%   verbose  print the table (default true)
%
% OUTPUT struct array S, one entry per dlog variable:
%   .name, .R2 (share of variance explained by the month dummies),
%   .F (joint test statistic on the eleven dummies), .df ([11 T-12]),
%   .flagged (logical, R2 above the 0.10 rule of thumb)
%
% No p-value is reported so the function stays toolbox-free. As a rule of
% thumb, with ~370 monthly observations an F above roughly 2.5 on (11, 358)
% degrees of freedom is well past conventional significance.

    if nargin < 2 || isempty(verbose), verbose = true; end

    Y = tr.Y;  dts = tr.dates;  names = tr.names;  tcode = tr.tcode;
    [T, K] = size(Y);
    dv = datevec(dts);  mo = dv(:,2);

    % Constant plus dummies for February..December.
    X = ones(T, 12);
    for m = 2:12
        X(:, m) = double(mo == m);
    end
    dfNum = 11;  dfDen = T - 12;

    S = struct('name',{},'R2',{},'F',{},'df',{},'flagged',{});
    for j = 1:K
        if ~strcmp(tcode{j}, 'dlog'), continue; end
        y  = Y(:,j);
        b  = X \ y;
        e  = y - X*b;
        yc = y - mean(y);
        R2 = 1 - (e.'*e) / (yc.'*yc);
        F  = (R2/dfNum) / ((1-R2)/dfDen);
        n = numel(S) + 1;
        S(n).name    = names{j};
        S(n).R2      = R2;
        S(n).F       = F;
        S(n).df      = [dfNum dfDen];
        S(n).flagged = R2 > 0.10;
    end

    if ~verbose, return; end

    fprintf('\n--- Deterministic monthly seasonality in the growth rates ---\n');
    fprintf('  %-8s %8s %10s\n', 'series', 'R2', sprintf('F(11,%d)', dfDen));
    for n = 1:numel(S)
        mark = '';
        if S(n).flagged, mark = '  <-- seasonal'; end
        fprintf('  %-8s %8.3f %10.1f%s\n', S(n).name, S(n).R2, S(n).F, mark);
    end
    if any([S.flagged])
        fprintf(2, ['  Seasonal variance is orthogonal to the rest of the system, so it\n' ...
                    '  inflates that variable''s own share in the decomposition and shrinks\n' ...
                    '  every other shock''s. Substitute a seasonally adjusted series (for\n' ...
                    '  CPI: opts.cpiFile) or apply X-13 to the level before differencing.\n']);
    end
    fprintf('\n');
end
