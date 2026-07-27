%RUN_SOFT_CONDITIONS  The "Soft conditions" segment of Running_code.txt.
%
% Waggoner-Zha soft conditioning: instead of pinning the constrained variable to
% an exact path, condition on a DISTRIBUTION around it by adding a .sd field to
% any spec. That gives the constrained variable its own band -- more honest for
% something labelled a risk scenario. .sd may be a scalar, or a per-horizon
% vector (tight early, loose later).
%   * Tightening, roughly            (policy path, sd 0.4)
%   * Oil spike, loosely held        (pulse, sd 0.03 in dlog units)
%   * Unemp ~5%, 1yr (soft)          (hold + window + sd 0.5)
%   * Rate +100bp, tight then loose  (per-horizon sd vector)
%   * CAD ~1.42 (soft)               (ramp, sd 0.012 in dlog units)
%
% Figures go to ./figures/soft_conditions.
%
% USAGE (from the project folder)
%   >> run_soft_conditions
% Rs is left in the workspace.

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

os = ob;
os.figSub = 'soft_conditions';
os.scenarios = { ...
  struct('var','rate','type','policy','step',0.25,'peak',4,'hold',8,'floor',2.5, ...
         'sd',0.4, 'name','Tightening, roughly'), ...
  struct('var','oil','type','pulse','to',110,'months',3,'hold',3, ...
         'sd',0.03,'name','Oil spike, loosely held'), ...
  struct('var','unemp','type','hold','at',5.0,'window',[1 12],'sd',0.5, ...
         'name','Unemp ~5%, 1yr (soft)'),                                    ... % soft + window
  struct('var','rate','type','ramp','to',L(ix.rate)+1.0,'months',6, ...
         'sd',[0.1*ones(6,1); 0.5*ones(ob.H-6,1)], ...
         'name','Rate +100bp, tight then loose'),                           ... % per-horizon sd
  struct('var','fx','type','ramp','to',1.42,'months',4,'sd',0.012, ...
         'name','CAD ~1.42 (soft)')                                         ... % soft, dlog units
};

Rs = main_risk_scenarios('fred', os);
fprintf('run_soft_conditions: done. Rs is in the workspace.\n');
