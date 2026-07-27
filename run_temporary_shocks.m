%RUN_TEMPORARY_SHOCKS  The "Temporary Shocks" segment of Running_code.txt.
%
% Partial-horizon / transitory conditions: the constraint binds for a few months
% and then the variable is left FREE so the model, not the analyst, decides the
% tail. Uses the 'pulse' type and the .window field.
%   * Oil spike, 6 months            (pulse: ramp to $110 over 3m, hold 3m, free)
%   * Oil at $110, 6m then free      (hold pinned to months 1-6 only)
%   * Rate +150bp, 1yr then free     (temporary tightening, then normalises)
%   * Tight labour, 9m then free     (transitory labour tightness, horizons 1-9)
%   * CAD weakens, 1Q then free      (transitory depreciation, first quarter)
%   * Equity selloff, 6m then free   (short replay of a 2008 drawdown, then recovers)
%
% Figures go to ./figures/temporary_shocks.
%
% USAGE (from the project folder)
%   >> run_temporary_shocks
% Rp is left in the workspace.

%% --- put the project (root + src/) on the MATLAB path, from any working dir ---
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir); addpath(fullfile(thisDir,'src'));
assert(exist('main_risk_scenarios','file')==2, ...
    'main_risk_scenarios.m not found under ./src -- keep this script at the project root.');

ob = scenario_config();

% origin levels (pass 1)
[L, ix, originMonth] = scenarioOrigin(ob);

ov = datevec(originMonth);
ob.figStart = [ov(1)-2, ov(2)];
ob.bands    = true;
ob.figures  = true;
ob.bandOpts = struct('nDraws', 400);

op = ob;
op.figSub = 'temporary_shocks';
op.scenarios = { ...
  % pulse: ramp to $110 over 3m, hold 3m, then FREE
  struct('var','oil','type','pulse','to',110,'months',3,'hold',3, ...
         'name','Oil spike, 6 months'), ...
  % .window: oil pinned at $110 for months 1-6 only, tail left to the model
  struct('var','oil','type','hold','at',110,'window',[1 6], ...
         'name','Oil at $110, 6m then free'), ...
  % temporary tightening that then normalises endogenously
  struct('var','rate','type','pulse','to',L(ix.rate)+1.5,'months',6,'hold',6, ...
         'name','Rate +150bp, 1yr then free'), ...
  % transitory labour tightness, horizons 1-9 only
  struct('var','unemp','type','hold','at',max(L(ix.unemp)-1.5,4.5),'window',[1 9], ...
         'name','Tight labour, 9m then free'), ...
  % transitory CAD depreciation, first quarter only (fx up = CAD weaker)
  struct('var','fx','type','ramp','to',1.45,'months',3,'window',[1 3], ...
         'name','CAD weakens, 1Q then free'), ...
  % replay a short equity drawdown, then let it recover freely
  struct('var','tsx','type','replay','from',[2008 9],'months',6, ...
         'name','Equity selloff, 6m then free') };

Rp = main_risk_scenarios('fred', op);
fprintf('run_temporary_shocks: done. Rp is in the workspace.\n');
