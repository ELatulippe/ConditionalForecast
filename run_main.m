%RUN_MAIN  Baseline run -- the "Main" segment of Running_code.txt.
%
% Estimates the monthly VAR on the compact panel (two extra variables: the BoC
% non-energy commodity index and the 10y GoC yield), draws the baseline fan
% charts, and returns everything in R.
%
% USAGE
%   Keep this script at the project root (it adds src/ to the path itself),
%   then from MATLAB/Octave:
%       >> run_main
%   R is left in the workspace.
%
% Requires network access (source = 'fred') and the StatCan CSVs in data/.

%% --- put the project (root + src/) on the MATLAB path, from any working dir ---
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir); addpath(fullfile(thisDir,'src'));
assert(exist('main_risk_scenarios','file')==2, ...
    'main_risk_scenarios.m not found under ./src -- keep this script at the project root.');

ob = struct();
ob.files     = struct('house','1810020501-eng.csv', ...
                      'cpi',  '1810000601-eng.csv', ...
                      'unemp','1410028701-eng.csv');
ob.rows      = struct('unemp','Unemployment rate');
ob.sample    = [1992 1];
ob.cacheFile = 'panel_k9.mat';        % different variable set -> its own cache
ob.covid     = [2020 4; 2020 6];
ob.seasonalDummies = true;            % bcne is seasonal (R2 = 0.12, F = 5.0)
ob.gate      = false;
ob.bands     = true;
ob.bandOpts  = struct('nDraws', 400);
ob.H         = 48;                    % estimate over four years
ob.figStart  = [2024 1];
ob.figEnd    = [2028 6];              % but plot two
ob.paperFigs = false;
ob.diffFigs  = false;
ob.p         = 12;
ob.prior     = struct('lambda', 0.2);

ob.extra = { ...
  struct('name','bcne','file','BCPI_MONTHLY-sd-1972-01-01.csv', ...
         'row','M.BCNE','tcode','dlog','after','oil'), ...
  struct('name','long','id','IRLTLT01CAM156N', ...
         'tcode','level','after','rate') };

R = main_risk_scenarios('fred', ob);
fprintf('run_main: done. Figures in ./figures. Struct R is in the workspace.\n');
