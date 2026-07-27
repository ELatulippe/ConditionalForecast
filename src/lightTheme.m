function lightTheme(f)
%LIGHTTHEME  Force a white-background, dark-ink theme on a figure before export.
%
%   lightTheme(f)     % f a figure handle (default gcf)
%
% Call this just before print/saveas so the saved PNG is light (for the paper
% and the HTML report) even if a dark default theme is in effect. It overrides,
% for figure f only, whatever background/foreground colours are currently set:
% the figure and every axes background become white, axis lines and tick labels
% become dark, the grid becomes light grey, any pure-white text (dark-theme
% labels/titles) is recoloured dark, and InvertHardcopy is left on so print()
% keeps the white background. Coloured series (the green history, red median,
% shaded bands) are untouched.
%
% Why this is the fix: a black figure comes from setting the figure/axes 'Color'
% to black (often with 'InvertHardcopy','off' so print does not re-whiten it).
% Setting those back to white here is all that is needed; doing it at save time
% means you do not have to hunt down where the dark theme was applied.

    if nargin < 1 || isempty(f), f = gcf; end
    set(f, 'Color', 'w', 'InvertHardcopy', 'on');

    ax = findall(f, 'type', 'axes');
    for i = 1:numel(ax)
        a = ax(i);
        set(a, 'Color', 'w');
        try, set(a, 'XColor', [0.15 0.15 0.15], 'YColor', [0.15 0.15 0.15]); catch, end
        try, set(a, 'GridColor', [0.85 0.85 0.85]); catch, end        % R2014b+
        t = get(a, 'Title');
        if ~isempty(t) && ishghandle(t), set(t, 'Color', [0.10 0.10 0.10]); end
    end

    % recolour any pure-white text (dark-theme titles/labels/legend text) to dark
    tx = findall(f, 'type', 'text');
    for i = 1:numel(tx)
        try
            if isequal(round(get(tx(i), 'Color')), [1 1 1])
                set(tx(i), 'Color', [0.10 0.10 0.10]);
            end
        catch
        end
    end

    lg = findall(f, 'type', 'legend');
    for i = 1:numel(lg)
        try, set(lg(i), 'Color', 'w', 'TextColor', [0.10 0.10 0.10], ...
                        'EdgeColor', [0.80 0.80 0.80]); catch, end
    end
end
