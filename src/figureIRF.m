function files = figureIRF(IRF, shocks, outdir, opts)
% FIGUREIRF  Plot Cholesky impulse responses, one figure per shock.
%
%   figureIRF(IRF)                              % every shock
%   figureIRF(IRF, {'rate','oil'}, 'figures')
%
% Each figure shows how every variable responds to one shock, with a zero
% line. Read the CPI panel for the price puzzle: under a contractionary
% policy shock the price level should fall, so a path sitting above zero over
% the first year or two is the puzzle, and the band tells you whether it is
% distinguishable from zero.
%
% Panel units follow the transform. Variables entering in log-differences
% were cumulated by irfCholesky, so they read as percent deviations of the
% LEVEL; level variables read in percentage points. The title records the
% normalisation, which is what makes the numbers comparable to a literature.
%
% INPUTS
%   IRF     output of irfCholesky
%   shocks  names or indices of the shocks to plot (default: all)
%   outdir  folder for the .png files (default 'figures')
%   opts    .color  RGB for the response line
%           .prefix filename prefix (default 'figureIRF')
%
% OUTPUT: cell array of the files written.

    if nargin < 2 || isempty(shocks), shocks = 1:numel(IRF.names); end
    if nargin < 3 || isempty(outdir), outdir = 'figures'; end
    if nargin < 4, opts = struct(); end
    if ~isfield(opts,'color'),  opts.color  = [0.20 0.35 0.70]; end
    if ~isfield(opts,'prefix'), opts.prefix = 'figureIRF'; end
    if ~exist(outdir,'dir'), mkdir(outdir); end

    if ischar(shocks), shocks = {shocks}; end
    if iscell(shocks)
        sel = zeros(1,numel(shocks));
        for i = 1:numel(shocks)
            j = find(strcmpi(shocks{i}, IRF.names), 1);
            if isempty(j)
                error('figureIRF:shock', 'No shock named "%s". Available: %s', ...
                      shocks{i}, strjoin(IRF.names, ', '));
            end
            sel(i) = j;
        end
    else
        sel = shocks;
    end

    K = numel(IRF.names);  H = IRF.horizon;
    nr = ceil(sqrt(K));  nc = ceil(K/nr);
    hasBand = isfield(IRF,'lo');
    files = {};

    for s = sel
        f = figure('Visible','off','Position',[100 100 1000 750]);
        for v = 1:K
            subplot(nr, nc, v); hold on; box on;
            y = squeeze(IRF.resp(v,s,:));
            if hasBand
                lo = squeeze(IRF.lo(v,s,:));  hi = squeeze(IRF.hi(v,s,:));
                fillBand((1:H).', lo, hi, opts.color, 0.22);
            end
            plot(1:H, y, '-', 'Color', opts.color, 'LineWidth', 1.6);
            plot([1 H], [0 0], 'k--', 'LineWidth', 0.9);
            title(sprintf('%s (%s)', IRF.names{v}, unitOf(IRF, v)));
            xlim([1 H]); grid on;
            if v > K - nc, xlabel('months'); end
        end
        annotateTitle(sprintf('Response to a %s shock%s', IRF.names{s}, normNote(IRF, s)));
        path = fullfile(outdir, sprintf('%s_%s.png', opts.prefix, ...
                        regexprep(lower(IRF.names{s}), '[^a-z0-9]+', '_')));
        saveFig(f, path);
        files{end+1} = path; %#ok<AGROW>
    end
end

% ======================================================================
function u = unitOf(IRF, v)
    if strcmp(IRF.tcode{v}, 'dlog')
        if IRF.cumulate, u = '% level'; else, u = '% growth'; end
    else
        u = 'pp';
    end
end

function t = normNote(IRF, s)
    if isempty(IRF.scaleTo)
        t = sprintf('  [1 sd = %.3f on %s]', IRF.scale(s), IRF.names{s});
    else
        t = sprintf('  [normalised to %+g on %s]', IRF.scaleTo, IRF.names{s});
    end
end

function annotateTitle(txt)
% A figure-level title without depending on sgtitle, which older MATLAB and
% Octave do not have.
    ax = axes('Position',[0 0.95 1 0.05],'Visible','off');
    text(0.5, 0.5, txt, 'Parent', ax, 'HorizontalAlignment','center', ...
         'FontWeight','bold', 'Interpreter','none');
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

function saveFig(f, path)
    lightTheme(f);            % force white background for the paper/report
    try
        print(f, path, '-dpng', '-r120');
    catch
        saveas(f, path);
    end
    close(f);
end
