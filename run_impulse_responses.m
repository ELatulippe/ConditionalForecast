%RUN_IMPULSE_RESPONSES  The "Impulse Responses" segment of Running_code.txt.
%
% Structural impulse responses under the recursive (Cholesky) ordering. This is
% a SELF-CONTAINED configuration: it uses a larger four-extra panel (bcne, long,
% emp, tsx) cached as panel_k13.mat, and its own p = 12 estimate with no prior
% and no scenario figures. irfCholesky rescales each shock so its own variable
% moves by 1 on impact (rate -> per 100bp), and figureIRF writes one figure per
% shock to ./figures/irf.
%
% USAGE (from the project folder)
%   >> run_impulse_responses
% R12 (the estimate) and I12 (the IRFs) are left in the workspace.

%% --- put the project (root + src/) on the MATLAB path, from any working dir ---
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir); addpath(fullfile(thisDir,'src'));
assert(exist('main_risk_scenarios','file')==2, ...
    'main_risk_scenarios.m not found under ./src -- keep this script at the project root.');
assert(exist('irfCholesky','file')==2, 'irfCholesky.m not found on the path.');

ob = struct();
ob.files     = struct('house','1810020501-eng.csv', ...
                      'cpi',  '1810000601-eng.csv', ...
                      'unemp','1410028701-eng.csv');
ob.rows      = struct('unemp','Unemployment rate');
ob.sample    = [1992 1];
ob.cacheFile = 'panel_k13.mat';           % was panel_k9 -- bumped for the larger panel
ob.covid     = [2020 4; 2020 6];
ob.seasonalDummies = true;
ob.gate      = false;
ob.figures   = false;
ob.p         = 12;
ob.extra = { ...
  struct('name','bcne','file','BCPI_MONTHLY-sd-1972-01-01.csv', ...
         'row','M.BCNE','tcode','dlog','after','oil'), ...
  struct('name','long', 'id','IRLTLT01CAM156N','tcode','level','after','rate'),  ... % 10y GoC
  struct('name','emp',  'id','LFEMTTTTCAM647S', 'tcode','dlog', 'after','unemp'), ... % CA employment
  struct('name','tsx',  'id','SPASTT01CAM661N', 'tcode','dlog', 'after','long')   ... % equity, late
};

R12 = main_risk_scenarios('fred', ob);
I12 = irfCholesky(R12, 48, struct('scaleTo', 1, 'boot', 400));
figureIRF(I12, {'rate','oil'}, fullfile('figures','irf'), struct('prefix','irf_k13_p12'));
fprintf('run_impulse_responses: done. Figures in ./figures/irf.\n');
