%RUN_VARIANCE_DECOMP  The "Variance Decomposition" segment of Running_code.txt.
%
% Forecast-error variance decomposition of the estimated VAR. Two outputs:
%   (1) a LaTeX fragment (figures/fevd_rows.tex) with the share of each target's
%       FEV due to the oil / usip / unemp / rate shocks at horizons 3,12,24,48;
%   (2) grouped stacked-area charts (figures/fevd/fevd_<target>.png) collapsing
%       the shocks into Commodity / US / CA-prices / CA-labour / CA-finpol blocks.
%
% USAGE (from the project folder)
%   >> run_variance_decomp
%
% Estimates the VAR itself (via scenarioOrigin), so it is independent of the
% other run_*.m scripts. "Rb" below is that estimate.

%% --- put the project (root + src/) on the MATLAB path, from any working dir ---
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir); addpath(fullfile(thisDir,'src'));
assert(exist('main_risk_scenarios','file')==2, ...
    'main_risk_scenarios.m not found under ./src -- keep this script at the project root.');
assert(exist('varianceDecomp','file')==2, 'varianceDecomp.m not found on the path.');

% Estimate once to get the model + variable index (the note's "Rb").
ob = scenario_config();
[~, ~, ~, Rb] = scenarioOrigin(ob);

[shares, MSPE] = varianceDecomp(Rb.model, 48, [], [3 12 24 48], Rb.data.names, Rb.data.names);
ix = Rb.data.idx;
shocks  = {'oil','usip','unemp','rate'};
targets = {'cpi','house','unemp','rate','fx'};
Hrep    = [3 12 24 48];

%% ---- LaTeX table fragment ----
if ~exist('figures','dir'), mkdir('figures'); end
fid = fopen(fullfile('figures','fevd_rows.tex'),'w');
for s = 1:numel(shocks)
  fprintf(fid,'\\multicolumn{%d}{@{}l}{\\textit{Panel %c: %s shock}}\\\\\n', ...
          numel(targets)+1, 'A'+s-1, shocks{s});
  for hh = Hrep
    v = 100*arrayfun(@(t) shares(ix.(targets{t}), ix.(shocks{s}), hh), 1:numel(targets));
    fprintf(fid,'%d', hh); fprintf(fid,' & %.0f', v); fprintf(fid,'\\\\\n');
  end
  fprintf(fid,'\\addlinespace\n');
end
fclose(fid);

%% ---- grouped stacked-area charts ----
Hn = 48;
% Each shock group, with the members that make it up and a CLEAR legend label
% (so "US" etc. is unambiguous). members must be variable names present in ix.
groups = { ...
  struct('members',{{'oil','bcne'}},          'label','Commodity (WTI, non-energy)'), ...
  struct('members',{{'usip','uscpi','usrate'}},'label','United States (IP, CPI, T-bill)'), ...
  struct('members',{{'cpi','house'}},          'label','Canada prices (CPI, housing)'), ...
  struct('members',{{'unemp','emp'}},          'label','Canada labour (unemp, employment)'), ...
  struct('members',{{'rate','long','tsx','fx'}},'label','Canada finance/policy (rate, 10y, TSX, FX)') };
nG      = numel(groups);
% Keep only groups whose members are actually in the panel (robust to spec
% changes), so the legend never lists a group with nothing behind it.
keep = false(1,nG);
for g = 1:nG
  keep(g) = all(cellfun(@(nm) isfield(ix,nm), groups{g}.members));
end
groups = groups(keep);  nG = numel(groups);
glabels = cellfun(@(gg) gg.label, groups, 'UniformOutput', false);
outdir = fullfile('figures','fevd'); if ~exist(outdir,'dir'), mkdir(outdir); end

for tv = {'cpi','house','unemp','rate','fx'}
  k = ix.(tv{1});
  G = zeros(Hn, nG);
  for g = 1:nG
    cols = cellfun(@(nm) ix.(nm), groups{g}.members);
    G(:,g) = 100*squeeze(sum(shares(k, cols, 1:Hn), 2));
  end
  figure('Color','w','Position',[100 100 1000 480]);
  hAr = area(1:Hn, G, 'LineStyle','none');   % keep the handles so labels map 1:1
  grid on; ylim([0 100]); xlim([1 Hn]);
  set(gca,'FontSize',12,'LineWidth',1.1);
  xlabel('horizon (months)'); ylabel('share of forecast-error variance, %');
  title(sprintf('Variance decomposition of %s', Rb.data.names{k}));
  lg = legend(hAr, glabels, 'Location','eastoutside','FontSize',10, 'Box','on');
  try, title(lg, 'shock group'); catch, end   % legend title (drop on old Octave)
  lightTheme(gcf);
  print(gcf, fullfile(outdir, sprintf('fevd_%s.png', tv{1})), '-dpng','-r150');
end

fprintf('run_variance_decomp: done. Table in figures/fevd_rows.tex, charts in figures/fevd.\n');
