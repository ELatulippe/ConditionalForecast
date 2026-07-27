function out = build_report(figdir)
%BUILD_REPORT  Assemble every generated figure + the FEVD table into one report.
%
%   build_report                 % uses <root>/figures
%   build_report('/path/figs')   % or a folder you choose
%
% Produces, in order of what your setup supports:
%
%   1. report.html  -- ALWAYS. A single self-contained file with the images
%      embedded (base64) and the FEVD table rendered as an HTML table. It opens
%      in any browser and prints to PDF (Ctrl+P -> Save as PDF) with no extra
%      software. This is the reliable default and needs no LaTeX.
%   2. report.pdf   -- if your MATLAB is R2021b+ (uses exportgraphics), a PDF is
%      built directly from the figures, again with NO LaTeX.
%   3. report.tex   -- always written; compiled to PDF only if pdflatex is on
%      the PATH. This is the typeset version for those who have LaTeX.
%
% Sections follow the derivation note: baseline / single, joint, temporary,
% soft, counterfactuals, impulse responses, variance decomposition. Missing
% figures are skipped, so a partial run still produces a report.
%
% OUT is a struct with the paths actually written (.html, .pdf, .tex).

    if nargin < 1 || isempty(figdir)
        here = fileparts(mfilename('fullpath'));
        figdir = fullfile(here,'figures');
    end
    if ~exist(figdir,'dir')
        error('build_report:nofigs', 'No figures folder at %s. Run some exercises first.', figdir);
    end

    % Section title | subfolder under figures/ ('' = the figures/ root)
    sections = { ...
        'Baseline and single-variable scenarios', ''; ...
        'Joint (multi-variable) scenarios',        'joint_scenarios'; ...
        'Temporary / partial-horizon shocks',      'temporary_shocks'; ...
        'Soft (distributional) conditions',        'soft_conditions'; ...
        'Historical counterfactuals',              'counterfactual'; ...
        'Impulse responses',                       'irf'; ...
        'Forecast-error variance decomposition',   'fevd' };
    fevdTable = fullfile(figdir,'fevd_rows.tex');

    out = struct('html','', 'pdf','', 'tex','');

    % ---- 1. HTML (always) ----------------------------------------------
    out.html = writeHTML(figdir, sections, fevdTable);

    % ---- 2. native PDF via exportgraphics (no LaTeX) -------------------
    try
        out.pdf = nativePDF(figdir, sections, fevdTable);
    catch err
        fprintf(2, 'build_report: native PDF skipped (%s)\n', err.message);
        out.pdf = '';
    end

    % ---- 3. LaTeX source, compiled only if pdflatex exists -------------
    out.tex = writeTeX(figdir, sections, fevdTable);
    if isempty(out.pdf)
        p = compileTeX(figdir, out.tex);
        if ~isempty(p), out.pdf = p; end
    end

    % ---- final guidance -------------------------------------------------
    fprintf('\nbuild_report: wrote\n');
    if ~isempty(out.html), fprintf('  HTML : %s\n', out.html); end
    if ~isempty(out.tex),  fprintf('  TeX  : %s\n', out.tex);  end
    if ~isempty(out.pdf)
        fprintf('  PDF  : %s\n', out.pdf);
    else
        fprintf(['  PDF  : not produced (no exportgraphics and no pdflatex).\n' ...
                 '         Open report.html and use your browser''s Print -> Save as PDF.\n']);
    end
end

% ======================================================================
% HTML report (self-contained; the reliable, dependency-free output)
% ======================================================================
function htmlPath = writeHTML(figdir, sections, fevdTable)
    htmlPath = fullfile(figdir,'report.html');
    fid = fopen(htmlPath,'w');
    if fid < 0, error('build_report:open','Cannot write %s', htmlPath); end
    oc = onCleanup(@() fclose(fid));

    fprintf(fid, '<!doctype html><html><head><meta charset="utf-8">\n');
    fprintf(fid, '<title>Canadian Risk Scenarios - Results</title>\n');
    fprintf(fid, ['<style>body{font-family:Segoe UI,Arial,sans-serif;max-width:960px;' ...
                  'margin:24px auto;color:#1a1a1a;padding:0 16px}h1{font-size:24px}' ...
                  'h2{font-size:18px;border-bottom:2px solid #ccc;padding-bottom:4px;' ...
                  'margin-top:32px}figure{margin:18px 0;text-align:center}' ...
                  'img{max-width:100%%;height:auto;border:1px solid #e2e2e2}' ...
                  'figcaption{font-size:13px;color:#555;margin-top:6px}' ...
                  'table{border-collapse:collapse;font-size:13px;margin:8px 0}' ...
                  'td,th{border:1px solid #ccc;padding:3px 8px;text-align:right}' ...
                  'td.panel{text-align:left;font-style:italic;background:#f4f4f4}' ...
                  '@media print{h2{page-break-before:always}}</style>\n']);
    fprintf(fid, '</head><body>\n');
    fprintf(fid, '<h1>Canadian Risk Scenarios &mdash; Results</h1>\n');
    fprintf(fid, '<p>Generated %s. Conditional-forecast (Waggoner&ndash;Zha) replication.</p>\n', ...
            htmlEscape(datestr(now,'yyyy-mm-dd')));

    for s = 1:size(sections,1)
        title = sections{s,1};  sub = sections{s,2};
        d = figdir; if ~isempty(sub), d = fullfile(figdir, sub); end
        pngs = listPNGs(d);
        isFEVD = strcmp(sub,'fevd');
        if isempty(pngs) && ~(isFEVD && exist(fevdTable,'file')), continue; end

        fprintf(fid, '<h2>%s</h2>\n', htmlEscape(title));
        if isFEVD && exist(fevdTable,'file')
            writeFEVDhtml(fid, fevdTable);
        end
        for i = 1:numel(pngs)
            uri = dataURI(pngs{i});
            fprintf(fid, '<figure><img src="%s"><figcaption>%s</figcaption></figure>\n', ...
                    uri, htmlEscape(humanize(pngs{i})));
        end
    end
    fprintf(fid, '</body></html>\n');
end

function writeFEVDhtml(fid, fevdTable)
    fprintf(fid, ['<table><tr><th>h</th><th>cpi</th><th>house</th><th>unemp</th>' ...
                  '<th>rate</th><th>fx</th></tr>\n']);
    frag = fopen(fevdTable,'r');
    if frag < 0, fprintf(fid,'</table>\n'); return; end
    ln = fgetl(frag);
    while ischar(ln)
        t = strtrim(ln);
        if isempty(t) || ~isempty(strfind(t,'addlinespace'))          %#ok<STREMP>
            ln = fgetl(frag); continue;
        end
        t = strrep(t, '\\', '');                                       % drop trailing row break
        if ~isempty(strfind(t,'multicolumn'))                          %#ok<STREMP>
            tok = regexp(t, '\\textit\{([^}]*)\}', 'tokens', 'once');
            if isempty(tok), tok = {'Panel'}; end
            fprintf(fid, '<tr><td class="panel" colspan="6">%s</td></tr>\n', htmlEscape(tok{1}));
        elseif ~isempty(strfind(t,'&'))                                %#ok<STREMP>
            parts = strtrim(strsplit(t, '&'));
            fprintf(fid, '<tr>');
            for j = 1:numel(parts), fprintf(fid, '<td>%s</td>', htmlEscape(parts{j})); end
            fprintf(fid, '</tr>\n');
        end
        ln = fgetl(frag);
    end
    fclose(frag);
    fprintf(fid, '</table>\n<p style="font-size:12px;color:#666">Share of forecast-error variance (%%) by shock, horizons 3/12/24/48 months.</p>\n');
end

function uri = dataURI(pngPath)
% Embed the PNG as a base64 data URI so report.html is self-contained. Falls
% back to a relative path if base64 encoding is unavailable.
    try
        fid = fopen(pngPath,'r');
        raw = fread(fid, Inf, '*uint8'); fclose(fid);
        b64 = matlab.net.base64encode(raw);
        uri = ['data:image/png;base64,' b64];
    catch
        [~,nm,ext] = fileparts(pngPath);
        pd = fileparts(pngPath); [~,sub] = fileparts(pd);
        uri = [sub '/' nm ext];                    % relative to figures/
    end
end

% ======================================================================
% Native PDF via exportgraphics (R2021b+), no LaTeX, no toolboxes.
% ======================================================================
function pdfPath = nativePDF(figdir, sections, fevdTable)
    pdfPath = '';
    if exist('exportgraphics','file') ~= 2, return; end
    if exist('verLessThan','file') == 2 && verLessThan('matlab','9.11'), return; end

    pdfPath = fullfile(figdir,'report.pdf');
    if exist(pdfPath,'file'), delete(pdfPath); end
    f = figure('Color','w','Visible','off','Units','pixels','Position',[100 100 1000 720]);
    oc = onCleanup(@() close(f));

    % title page
    textPage(f, {'Canadian Risk Scenarios', 'Results', '', datestr(now,'yyyy-mm-dd')}, 20);
    exportgraphics(f, pdfPath, 'Append', true, 'Resolution', 150);

    for s = 1:size(sections,1)
        title = sections{s,1};  sub = sections{s,2};
        d = figdir; if ~isempty(sub), d = fullfile(figdir, sub); end
        pngs = listPNGs(d);
        isFEVD = strcmp(sub,'fevd');
        if isempty(pngs) && ~(isFEVD && exist(fevdTable,'file')), continue; end

        textPage(f, {title}, 18);
        exportgraphics(f, pdfPath, 'Append', true, 'Resolution', 150);

        if isFEVD && exist(fevdTable,'file')
            textPage(f, fevdLines(fevdTable), 10, 'FixedWidth');
            exportgraphics(f, pdfPath, 'Append', true, 'Resolution', 150);
        end
        for i = 1:numel(pngs)
            clf(f);
            ax = axes(f, 'Position',[0.03 0.05 0.94 0.88]);
            image(ax, imread(pngs{i})); axis(ax, 'image', 'off');
            title(ax, humanize(pngs{i}), 'Interpreter','none', 'FontSize',12);
            exportgraphics(f, pdfPath, 'Append', true, 'Resolution', 150);
        end
    end
end

function textPage(f, lines, fs, fontname)
    if nargin < 4, fontname = 'Helvetica'; end
    clf(f);
    ax = axes(f, 'Position',[0.06 0.05 0.88 0.9]); axis(ax,'off');
    y = 0.95;
    for i = 1:numel(lines)
        text(ax, 0.0, y, lines{i}, 'FontSize', fs, 'FontName', fontname, ...
             'Interpreter','none', 'VerticalAlignment','top');
        y = y - 0.035;
    end
end

function L = fevdLines(fevdTable)
    L = {'Variance decomposition (share of FEV, %), horizons 3/12/24/48', ...
         'h    cpi  house  unemp  rate   fx', ''};
    frag = fopen(fevdTable,'r'); if frag < 0, return; end
    ln = fgetl(frag);
    while ischar(ln)
        t = strtrim(ln);
        if ~isempty(t) && isempty(strfind(t,'addlinespace'))           %#ok<STREMP>
            t = strrep(t,'\\','');
            if ~isempty(strfind(t,'multicolumn'))                      %#ok<STREMP>
                tok = regexp(t, '\\textit\{([^}]*)\}', 'tokens', 'once');
                if ~isempty(tok), L{end+1} = tok{1}; end                %#ok<AGROW>
            elseif ~isempty(strfind(t,'&'))                            %#ok<STREMP>
                p = strtrim(strsplit(t,'&'));
                L{end+1} = sprintf('%-4s %4s %6s %6s %5s %5s', p{:}); %#ok<AGROW>
            end
        end
        ln = fgetl(frag);
    end
    fclose(frag);
end

% ======================================================================
% LaTeX source (compiled only if pdflatex is present)
% ======================================================================
function texPath = writeTeX(figdir, sections, fevdTable)
    texPath = fullfile(figdir,'report.tex');
    fid = fopen(texPath,'w');
    if fid < 0, error('build_report:open','Cannot write %s', texPath); end
    oc = onCleanup(@() fclose(fid));

    fprintf(fid, '%% Auto-generated by build_report.m -- do not edit by hand.\n');
    fprintf(fid, '\\documentclass[11pt]{article}\n');
    fprintf(fid, '\\usepackage[margin=1in]{geometry}\n');
    fprintf(fid, '\\usepackage{graphicx,booktabs,float,caption}\n');
    fprintf(fid, '\\captionsetup{font=small}\n');
    fprintf(fid, '\\title{Canadian Risk Scenarios --- Results}\n');
    fprintf(fid, '\\author{Conditional Forecast (Waggoner--Zha) replication}\n');
    fprintf(fid, '\\date{%s}\n', texEscape(datestr(now,'yyyy-mm-dd')));
    fprintf(fid, '\\begin{document}\n\\maketitle\n');

    for s = 1:size(sections,1)
        title = sections{s,1};  sub = sections{s,2};
        d = figdir; if ~isempty(sub), d = fullfile(figdir, sub); end
        pngs = listPNGs(d);
        isFEVD = strcmp(sub,'fevd');
        if isempty(pngs) && ~(isFEVD && exist(fevdTable,'file')), continue; end

        fprintf(fid, '\\section*{%s}\n', texEscape(title));
        % Inline the FEVD fragment (do NOT \input it inside a tabular -- that
        % disturbs the alignment and makes the leading \multicolumn "misplaced").
        if isFEVD && exist(fevdTable,'file')
            fprintf(fid, '\\begin{table}[H]\\centering\\small\n');
            fprintf(fid, '\\caption{Share of forecast-error variance (\\%%) by shock, at horizons 3/12/24/48 months.}\n');
            fprintf(fid, '\\begin{tabular}{l rrrrr}\n\\toprule\n');
            fprintf(fid, 'h & cpi & house & unemp & rate & fx \\\\\n\\midrule\n');
            frag = fopen(fevdTable,'r');
            if frag >= 0
                ln = fgetl(frag);
                while ischar(ln), fprintf(fid, '%s\n', ln); ln = fgetl(frag); end
                fclose(frag);
            end
            fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n');
        end
        for i = 1:numel(pngs)
            rel = relPath(figdir, pngs{i});
            fprintf(fid, '\\begin{figure}[H]\\centering\n');
            fprintf(fid, '\\includegraphics[width=0.92\\linewidth]{\\detokenize{%s}}\n', rel);
            fprintf(fid, '\\caption{%s}\n', texEscape(humanize(pngs{i})));
            fprintf(fid, '\\end{figure}\n');
        end
        fprintf(fid, '\\clearpage\n');
    end
    fprintf(fid, '\\end{document}\n');
end

function pdfPath = compileTeX(figdir, texPath)
    pdfPath = '';
    [st,~] = system('pdflatex --version');
    if st ~= 0, return; end
    [~, base] = fileparts(texPath);
    cmd = sprintf('cd "%s" && pdflatex -interaction=nonstopmode -halt-on-error "%s" > report.log 2>&1', ...
                  figdir, [base '.tex']);
    system(cmd);  system(cmd);
    cand = fullfile(figdir, [base '.pdf']);
    if exist(cand,'file'), pdfPath = cand; end
end

% ======================================================================
% shared helpers
% ======================================================================
function pngs = listPNGs(d)
    pngs = {};
    if ~exist(d,'dir'), return; end
    L = dir(fullfile(d,'*.png'));
    L = L(~[L.isdir]);
    [~,order] = sort({L.name});
    for i = order(:)'
        pngs{end+1} = fullfile(d, L(i).name); %#ok<AGROW>
    end
end

function r = relPath(base, f)
    base = strrep(base, '\', '/');  f = strrep(f, '\', '/');
    if ~isempty(base) && base(end) ~= '/', base = [base '/']; end
    if strncmp(f, base, numel(base)), r = f(numel(base)+1:end); else, r = f; end
end

function s = humanize(pngPath)
    [~, name] = fileparts(pngPath);
    for pre = {'figure5_bands_','figure_','fevd_','irf_'}
        if strncmp(name, pre{1}, numel(pre{1})), name = name(numel(pre{1})+1:end); end
    end
    s = strtrim(strrep(name, '_', ' '));
    if ~isempty(s), s(1) = upper(s(1)); end
end

function s = texEscape(s)
    s = strrep(s, '\', '\textbackslash{}');
    for c = {'&','%','$','#','_','{','}'}
        s = strrep(s, c{1}, ['\' c{1}]);
    end
    s = strrep(s, '~', '\textasciitilde{}');
    s = strrep(s, '^', '\textasciicircum{}');
end

function s = htmlEscape(s)
    s = strrep(s, '&', '&amp;');
    s = strrep(s, '<', '&lt;');
    s = strrep(s, '>', '&gt;');
    s = strrep(s, '"', '&quot;');
end
