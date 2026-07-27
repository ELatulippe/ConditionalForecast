function figureBands(data, idx, Fbase, B, k, H, outdir, xlo, showBands, xhi)
% FIGUREBANDS  Scenario forecast figure, median lines only by default.
%
%   figureBands(data, idx, Fbase, B, k, H, outdir)
%   figureBands(data, idx, Fbase, B, k, H, outdir, xlo)
%   figureBands(data, idx, Fbase, B, k, H, outdir, xlo, false)  % hide bands
%
% Draws the five Canadian variables with history, two forecast lines and the
% 68% / 90% bands around the scenario.
%
% ONE FORECAST LINE: the median of the PREDICTIVE distribution under the
% scenario, across bootstrap draws.
%
% No baseline is drawn. Neither candidate belongs on these axes. The OLS point
% baseline is not the centre of the bootstrap distribution -- the DGP is
% bias-adjusted and explosive replications are rejected, which shifts the
% baseline substantially for the persistent variables and negligibly for the
% transitory ones -- so plotting it next to a bootstrap median invites reading
% that shift as a scenario effect. A bootstrap median baseline avoids that but
% still invites comparing the gap by eye, which is a difference of medians
% rather than the median of differences. figureDiffBands does the comparison
% properly, within each draw.
%
% The shaded band is PREDICTIVE: where the variable might land under the
% scenario, future shocks included. It is not a test of whether the scenario
% differs from the baseline -- most of its width is shock uncertainty that the
% baseline shares. The vertical gap between the two lines is close to the
% scenario effect but not exactly it, since a median of a sum is not the sum
% of medians under the nonlinear display transforms. For the exact object, and
% for a band that answers the distinguishability question, use
% figureDiffBands, which differences within each draw.
%
% Every draw is transformed with dispSeries BEFORE the median is taken. The
% level and YoY transforms are nonlinear, so a median computed in the VAR's
% units and then transformed would not be the median of the displayed series.
%
% The constrained variable's panel shows all three lines collapsing onto the
% imposed path, which is correct: a hard condition fixes it exactly.
%
% INPUTS
%   data, idx   from loadCanadaData
%   Fbase       H x K baseline point forecast; used only to build the display
%               grid and the history segment, never plotted as a forecast
%   B           output of wzBands, with .condDraws and .baseMean present
%   k           which scenario to draw
%   H           horizon
%   outdir      folder for the .png
%   xlo         optional left-hand x limit (datenum)
%   showBands   optional, default true

    if ~isfield(B, 'condDraws')
        error('figureBands:noDraws', ...
            ['B has no .condDraws. Re-run wzBands with keepDraws = true; the ' ...
             'stored quantiles are in the VAR''s units and cannot be ' ...
             'transformed after the fact.']);
    end
    if nargin < 8, xlo = []; end
    if nargin < 9 || isempty(showBands), showBands = true; end

    specs = { idx.cpi,  'yoy',   'Prices (YoY %)';
              idx.house,'yoy',   'Housing Prices (YoY %)';
              idx.unemp,'level', 'Unemployment (%)';
              idx.rate, 'level', 'BC Rate (%)';
              idx.fx,   'yoy',   'Exchange Rate (YoY %)' };

    f = figure('Visible','off','Position',[100 100 1000 750]);
    for i = 1:5
        subplot(3,2,i);
        drawMedian(data, specs{i,1}, specs{i,2}, Fbase, B, k, H, specs{i,3}, ...
                   xlo, showBands, xhi);
    end

    subplot(3,2,6); axis off; hold on;
    col = [0.75 0.20 0.20];
    plot(nan,nan,'-','Color',[0.15 0.6 0.25],'LineWidth',1.5);
    plot(nan,nan,'-','Color',col,'LineWidth',1.6);
    lab = {'History','Scenario median'};
    if showBands
        fill([nan nan],[nan nan],col,'FaceAlpha',0.30,'EdgeColor','none');
        fill([nan nan],[nan nan],col,'FaceAlpha',0.15,'EdgeColor','none');
        lab = [lab, {'68% band','90% band'}];
    end
    legend(lab, 'Location','west','Box','off');
    title(B.scen(k).name);

    nm = regexprep(lower(B.scen(k).name), '[^a-z0-9]+', '_');
    saveFig(f, fullfile(outdir, sprintf('figure5_bands_%s.png', nm)));
end

% ======================================================================
function drawMedian(data, v, mode, Fbase, B, k, H, ttl, xlo, showBands, xhi)
    last = data.dates(end);
    N    = size(B.condDraws, 4);

    [ax, ~] = dispSeries(data, v, mode, Fbase(:,v), H);
    fut = ax >= last;  nf = sum(fut);  xf = ax(fut);

    % Transform each draw first, then take the median (and bands if asked).
    Cm = zeros(nf, N);
    for d = 1:N
        [~, c] = dispSeries(data, v, mode, B.condDraws(:,v,k,d), H);
        Cm(:,d) = c(fut);
    end
    qs = [0.05 0.16 0.50 0.84 0.95];
    Q  = zeros(nf, numel(qs));
    for t = 1:nf, Q(t,:) = qtile(Cm(t,:), qs); end

    [bx, bval] = dispSeries(data, v, mode, Fbase(:,v), H);
    hist = bx <= last;
    col  = [0.75 0.20 0.20];

    hold on; box on;
    if showBands
        fillBand(xf, Q(:,1), Q(:,5), col, 0.15);
        fillBand(xf, Q(:,2), Q(:,4), col, 0.30);
    end
    plot(bx(hist), bval(hist), '-', 'Color',[0.15 0.6 0.25], 'LineWidth',1.3);
    plot(xf, Q(:,3), '-', 'Color', col, 'LineWidth',1.6);

    yl = ylim; plot([last last], yl, '-', 'Color',[0.6 0.6 0.6]);
    title(ttl); grid on;
    if isempty(xlo), xlo = addMonthsLocal(last,-47); end
    if nargin < 10 || isempty(xhi), xhi = addMonthsLocal(last,H); end
    xlim([xlo, xhi]);
    datetick('x','yyyy','keeplimits');
end

% ======================================================================
function fillBand(x, lo, hi, col, alpha)
    ok = ~isnan(lo) & ~isnan(hi);
    if ~any(ok), return; end
    x = x(ok); lo = lo(ok); hi = hi(ok);
    try
        fill([x(:); flipud(x(:))], [lo(:); flipud(hi(:))], col, ...
             'FaceAlpha', alpha, 'EdgeColor','none');
    catch
        plot(x, lo, '--', 'Color', col);  plot(x, hi, '--', 'Color', col);
    end
end

function dn = addMonthsLocal(d, k)
    dv = datevec(d);
    n  = dv(1)*12 + (dv(2)-1) + k;
    dn = datenum(floor(n/12), mod(n,12)+1, 1);
end

function q = qtile(x, ps)
    x = sort(x(:));  n = numel(x);
    if n == 1, q = repmat(x, 1, numel(ps)); return; end
    grid = ((1:n) - 0.5) / n;
    q = zeros(1, numel(ps));
    for j = 1:numel(ps)
        pj = ps(j);
        if pj <= grid(1),       q(j) = x(1);
        elseif pj >= grid(end), q(j) = x(end);
        else
            i = find(grid <= pj, 1, 'last');
            w = (pj - grid(i)) / (grid(i+1) - grid(i));
            q(j) = (1-w)*x(i) + w*x(i+1);
        end
    end
end

function saveFig(f, path)
    lightTheme(f);            % force white background for the paper/report
    try
        print(f, path, '-dpng', '-r120');
    catch
        saveas(f, path);
    end
    close(f);
end
