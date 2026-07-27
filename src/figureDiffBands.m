function figureDiffBands(data, idx, B, k, H, outdir, xlo, xhi)
% FIGUREDIFFBANDS  Scenario minus baseline, with an uncertainty band.
%
%   figureDiffBands(data, idx, B, k, H, outdir)
%   figureDiffBands(data, idx, B, k, H, outdir, xlo)
%
% This is the figure that answers "is the scenario distinguishable from the
% baseline". Where the band excludes zero, the scenario moves the variable by
% more than the estimation uncertainty; where it straddles zero, it does not.
%
% WHY NOT JUST EYEBALL THE FAN CHART. In figureBands the black line is the
% baseline computed once from the OLS estimate, while the red line is the
% median across bootstrap draws. Those are not like for like: the bootstrap
% DGP is bias-adjusted and explosive replications are rejected, so the
% baseline itself shifts across draws and its bootstrap median sits away from
% the OLS point forecast. A visible gap between black and red therefore mixes
% the scenario effect with that shift, and can appear even for a scenario that
% does nothing.
%
% Here the difference is taken WITHIN each draw, from .condMean and .baseMean,
% so whatever moves the baseline moves the conditional forecast with it and
% cancels. What remains is the scenario effect and the parameter uncertainty
% around it. Shock uncertainty is absent by construction: with the parameters
% fixed the difference between the two conditional means is deterministic.
%
% Each draw is transformed with dispSeries and differenced in DISPLAY units
% before quantiles are taken, since the level and YoY transforms are nonlinear.
%
% The constrained variable's panel shows the imposed deviation with a
% degenerate band, which is correct: a hard condition fixes that path.

    if ~isfield(B, 'condMean') || ~isfield(B, 'baseMean')
        error('figureDiffBands:noMeans', ...
            ['B has no .condMean/.baseMean. Re-run wzBands with ' ...
             'keepDraws = true.']);
    end
    if nargin < 7, xlo = []; end

    specs = { idx.cpi,  'yoy',   'Prices (YoY %, diff)';
              idx.house,'yoy',   'Housing Prices (YoY %, diff)';
              idx.unemp,'level', 'Unemployment (pp diff)';
              idx.rate, 'level', 'BC Rate (pp diff)';
              idx.fx,   'yoy',   'Exchange Rate (YoY %, diff)' };

    f = figure('Visible','off','Position',[100 100 1000 750]);
    for i = 1:5
        subplot(3,2,i);
        drawDiff(data, specs{i,1}, specs{i,2}, B, k, H, specs{i,3}, xlo, xhi);
    end
    subplot(3,2,6); axis off; hold on;
    plot(nan,nan,'-','Color',[0.20 0.35 0.70],'LineWidth',1.6);
    fill([nan nan],[nan nan],[0.20 0.35 0.70],'FaceAlpha',0.30,'EdgeColor','none');
    fill([nan nan],[nan nan],[0.20 0.35 0.70],'FaceAlpha',0.15,'EdgeColor','none');
    plot(nan,nan,'k--','LineWidth',1.0);
    legend({'Median difference','68% band','90% band','zero'}, ...
           'Location','west','Box','off');
    nm = regexprep(lower(B.scen(k).name), '[^a-z0-9]+', '_');
    saveFig(f, fullfile(outdir, sprintf('figure6_diff_%s.png', nm)));
end

% ======================================================================
function drawDiff(data, v, mode, B, k, H, ttl, xlo, xhi)
    last = data.dates(end);
    N    = size(B.condMean, 4);

    [ax, ~] = dispSeries(data, v, mode, B.baseMean(:,v,1), H);
    fut = ax >= last;  nf = sum(fut);
    Dm  = zeros(nf, N);
    for d = 1:N
        [~, c] = dispSeries(data, v, mode, B.condMean(:,v,k,d), H);
        [~, b] = dispSeries(data, v, mode, B.baseMean(:,v,d),   H);
        Dm(:,d) = c(fut) - b(fut);          % difference in DISPLAY units
    end

    qs = [0.05 0.16 0.50 0.84 0.95];
    Q  = zeros(nf, numel(qs));
    for t = 1:nf, Q(t,:) = qtile(Dm(t,:), qs); end

    xf  = ax(fut);
    col = [0.20 0.35 0.70];
    hold on; box on;
    fillBand(xf, Q(:,1), Q(:,5), col, 0.15);
    fillBand(xf, Q(:,2), Q(:,4), col, 0.30);
    plot(xf, Q(:,3), '-', 'Color', col, 'LineWidth',1.6);
    plot([xf(1) xf(end)], [0 0], 'k--', 'LineWidth',1.0);
    title(ttl); grid on;
    if isempty(xlo), xlo = last; end
    if nargin < 9 || isempty(xhi), xhi = addMonthsLocal(last,H); end
    xlim([max(xlo,last), xhi]);
    datetick('x','yyyy','keeplimits');
end

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
