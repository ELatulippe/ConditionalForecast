function ob = scenario_config()
%SCENARIO_CONFIG  Shared panel / estimation options for the risk-scenario suite.
%
%   ob = scenario_config();
%
% Returns the common options struct used by every scenario, counterfactual and
% FEVD script in this folder (single / joint / temporary / soft / counterfactual
% / variance-decomposition). It fixes the panel (the six extra variables), the
% estimation window, the lag order, the Minnesota prior and the COVID / seasonal
% handling -- i.e. everything that must be IDENTICAL across those exercises so
% they all share one estimated VAR.
%
% Scenario-specific fields (.scenarios, .bands, .figures, .figStart, .figSub)
% are left OFF here and set by each run_*.m script.
%
% This mirrors the "Single Scenarios" configuration block in Running_code.txt.

    ob = struct();

    ob.files = struct('house','1810020501-eng.csv', ...
                      'cpi',  '1810000601-eng.csv', ...
                      'unemp','1410028701-eng.csv');
    ob.rows  = struct('unemp','Unemployment rate');

    ob.sample    = [1992 1];        % [startY startM] -> run to the latest month
    ob.cacheFile = 'panel_k9.mat';  % assembled panel cache (delete to rebuild)
    ob.covid     = [2020 4; 2020 6];
    ob.seasonalDummies = true;      % bcne is seasonal (R2 = 0.12, F = 5.0)
    ob.gate      = false;           % draw figures even if the Table-1 gate fails
    ob.H         = 48;              % forecast horizon, months
    ob.figEnd    = [2028 6];        % x-axis end for the fan charts
    ob.paperFigs = false;
    ob.diffFigs  = false;
    ob.p         = 12;              % lag order
    ob.prior     = struct('lambda', 0.2);   % Minnesota prior

    % ---- the six extra variables, in estimation (recursive) order ----
    ob.extra = { ...
      struct('name','bcne',  'file','BCPI_MONTHLY-sd-1972-01-01.csv', ...
             'row','M.BCNE','tcode','dlog', 'after','oil'),   ... % BoC non-energy commodities
      struct('name','long',  'id','IRLTLT01CAM156N','tcode','level','after','rate'),  ... % 10y GoC
      struct('name','uscpi', 'id','CPIAUCSL',        'tcode','dlog', 'after','usip'),  ... % foreign block, near usip
      struct('name','usrate','id','TB3MS',           'tcode','level','after','uscpi'), ... % US 3m T-bill
      struct('name','emp',   'id','LFEMTTTTCAM647S', 'tcode','dlog', 'after','unemp'), ... % CA employment
      struct('name','tsx',   'id','SPASTT01CAM661N', 'tcode','dlog', 'after','long')   ... % equity, late
    };
end
