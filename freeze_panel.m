%FREEZE_PANEL  Build the offline panel caches with ONE online run.
%
%   >> freeze_panel
%
% Fetches every series from FRED (plus the StatCan/BoC CSVs in data/), assembles
% the two panels the exercises use, and writes them to cache/ as .mat files:
%
%     cache/panel_k9.mat    the 9-variable scenario / counterfactual / FEVD panel
%     cache/panel_k13.mat   the larger impulse-response panel
%
% After this has run once (with internet), every run_*.m script and the
% regression tests can be run with NO network access: main_risk_scenarios reads
% the matching cache and returns it without contacting FRED. To make that a hard
% guarantee -- an error rather than a silent refetch if a cache is missing --
% set  ob.offline = true  in the config (see README, "Offline mode").
%
% Re-run freeze_panel whenever you want to refresh the data to the latest month,
% or after changing a panel's variable set (the cache is keyed to the exact set
% of variables and is rebuilt automatically when it no longer matches).

%% --- put the project (root + src/) on the MATLAB path, from any working dir ---
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir); addpath(fullfile(thisDir,'src'));
assert(exist('main_risk_scenarios','file')==2, ...
    'main_risk_scenarios.m not found under ./src -- keep this script at the project root.');

fprintf('\n=== freeze_panel: building offline caches (needs internet) ===\n\n');

% ---- 1. the scenario / counterfactual / FEVD panel (panel_k9) ----
fprintf('[1/2] scenario panel -> cache/panel_k9.mat\n');
ob = scenario_config();
ob.bands = false; ob.figures = false;      % just assemble + estimate, no output
main_risk_scenarios('fred', ob);

% ---- 2. the larger impulse-response panel (panel_k13) ----
fprintf('\n[2/2] impulse-response panel -> cache/panel_k13.mat\n');
oi = struct();
oi.files     = struct('house','1810020501-eng.csv', ...
                      'cpi',  '1810000601-eng.csv', ...
                      'unemp','1410028701-eng.csv');
oi.rows      = struct('unemp','Unemployment rate');
oi.sample    = [1992 1];
oi.cacheFile = 'panel_k13.mat';
oi.covid     = [2020 4; 2020 6];
oi.seasonalDummies = true;
oi.gate      = false;
oi.figures   = false;
oi.p         = 12;
oi.extra = { ...
  struct('name','bcne','file','BCPI_MONTHLY-sd-1972-01-01.csv', ...
         'row','M.BCNE','tcode','dlog','after','oil'), ...
  struct('name','long', 'id','IRLTLT01CAM156N','tcode','level','after','rate'),  ...
  struct('name','emp',  'id','LFEMTTTTCAM647S', 'tcode','dlog', 'after','unemp'), ...
  struct('name','tsx',  'id','SPASTT01CAM661N', 'tcode','dlog', 'after','long')   ...
};
main_risk_scenarios('fred', oi);

fprintf('\n=== freeze_panel: done. Caches written to cache/. ===\n');
fprintf('You can now run everything offline (e.g. run_all, or the run_*.m scripts).\n');
