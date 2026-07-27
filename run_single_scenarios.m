%RUN_SINGLE_SCENARIOS  The "Single Scenarios" segment of Running_code.txt.
%
% Four one-variable risk scenarios calibrated to the current forecast origin:
%   * Oil +30%   (ramp)
%   * Tight labour market for a year   (unemployment hold)
%   * Gradual tightening, +150bp   (policy path)
%   * U.S. recession   (replay of 2008 U.S. industrial production)
%
% Two passes, exactly as in the note:
%   pass 1  estimate once with no bands / no figures, read the origin levels;
%   pass 2  build the scenarios off those levels and draw the fan charts.
%
% USAGE (from the project folder)
%   >> run_single_scenarios
% Leaves Rb (the scenario run) in the workspace -- the counterfactual and FEVD
% scripts can reuse the same estimated VAR.

%% --- put the project (root + src/) on the MATLAB path, from any working dir ---
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir); addpath(fullfile(thisDir,'src'));
assert(exist('main_risk_scenarios','file')==2, ...
    'main_risk_scenarios.m not found under ./src -- keep this script at the project root.');

ob = scenario_config();

% ---- pass 1: read the forecast origin, no bands, no figures ----
[L, ix, originMonth] = scenarioOrigin(ob);

% ---- pass 2: scenarios calibrated to that origin ----
ov = datevec(originMonth);
ob.figStart = [ov(1)-2, ov(2)];        % two years of history before the origin
ob.bands    = true;
ob.figures  = true;
ob.bandOpts = struct('nDraws', 400);

ob.scenarios = { ...
  struct('var','oil',  'type','ramp', 'to',round(1.30*L(ix.oil)), 'months',3, ...
         'name',sprintf('Oil +30%% (to $%d)', round(1.30*L(ix.oil)))), ...
  struct('var','unemp','type','hold', 'at',max(L(ix.unemp)-1.5, 4.5), 'window',[1 12], ...
         'name','Tight labour market, 1yr'), ...
  struct('var','rate', 'type','policy','step',0.25,'peak',L(ix.rate)+1.5,'hold',12, ...
         'floor',L(ix.rate), 'name','Gradual tightening (+150bp)'), ...
  struct('var','usip', 'type','replay','from',[2008 1], ...
         'name','U.S. recession') };

Rb = main_risk_scenarios('fred', ob);
fprintf('run_single_scenarios: done. Rb is in the workspace.\n');
