%RUN_COUNTERFACTUALS  The two in-sample counterfactuals of Running_code.txt.
%
% Historical (in-sample) Waggoner-Zha counterfactuals over 2022-23. Nothing in
% the WZ algebra restricts it to the forecast horizon: wzCounterfactual keeps
% every realised shock except the one you override.
%
%   (1) "What if the Bank had tightened?"  -- confine the change to the policy
%       shock and push the rate 100bp above what actually happened.
%   (2) "What if the 2022 commodity spike hadn't happened?" -- oil enters in
%       log-differences, so "hold the LEVEL at its Dec-2021 value" is a path of
%       ZEROS in transformed units (no monthly growth).
%
% Both now produce the SAME two-panel figure (price level + unemployment on
% subplots) written to ./figures/counterfactual.
%
% USAGE (from the project folder)
%   >> run_counterfactuals
%
% This script estimates the VAR itself (via scenarioOrigin), so it does not
% depend on any other run_*.m having been executed first.

%% --- put the project (root + src/) on the MATLAB path, from any working dir ---
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir); addpath(fullfile(thisDir,'src'));
assert(exist('main_risk_scenarios','file')==2, ...
    'main_risk_scenarios.m not found under ./src -- keep this script at the project root.');
assert(exist('wzCounterfactual','file')==2, 'wzCounterfactual.m not found on the path.');

% Estimate once to get the model + transformed data (this is the "Rb" the note
% refers to -- only .model and .tr are used below).
ob = scenario_config();
[~, ~, ~, Rb] = scenarioOrigin(ob);

Hcf = 24;
ixr = Rb.tr.idx.rate;  ixc = Rb.tr.idx.cpi;  ixu = Rb.tr.idx.unemp;
ixo = Rb.tr.idx.oil;
outdir = fullfile('figures','counterfactual');
if ~exist(outdir,'dir'), mkdir(outdir); end

%% =====================================================================
%  Counterfactual 1 -- "what if the Bank had tightened?"  (+100bp)
%  =====================================================================
t0 = find(Rb.tr.dates == datenum(2021,12,1));
assert(~isempty(t0), 'Dec-2021 not found in the sample -- cannot place the counterfactual origin.');
rActual = Rb.tr.Y(t0+1:t0+Hcf, ixr);

Ct = wzCounterfactual(Rb.model, Rb.tr.Y, Rb.tr.dates, [2021 12], ...
                      ixr, rActual + 1.0, Hcf, [], ixr);   % tighter BoC, policy shock only
d  = Ct.dates;

pA = 100*cumsum(Ct.actual(:,ixc));  pC = 100*cumsum(Ct.cf(:,ixc));   % price: dlog -> cumulate
uA = Ct.actual(:,ixu);              uC = Ct.cf(:,ixu);               % unemp: already a level

figure('Color','w','Position',[100 100 1200 500]);

subplot(1,2,1);                                    % ---- price level ----
fill([d; flipud(d)], [pA; flipud(pC)], [0.90 0.95 0.90], 'EdgeColor','none'); hold on  % drop on Octave
hA = plot(d, pA, '-',  'LineWidth', 3, 'Color', [0.15 0.15 0.15]);
hC = plot(d, pC, '--', 'LineWidth', 3, 'Color', [0.15 0.55 0.25]);
datetick('x','yyyy-mm','keeplimits'); xlim([d(1) d(end)]);
grid on; set(gca,'FontSize',12,'LineWidth',1.2);
ylabel('cumulative price level, %'); title('Price level');
legend([hA hC], {'actual','tighter BoC (+100bp)'}, 'Location','northwest','FontSize',11);
gP = pC(end)-pA(end);
text(d(end), pC(end), sprintf('  %+.2f pp', gP), 'FontSize',12, ...
     'Color',[0.15 0.55 0.25], 'VerticalAlignment','middle');

subplot(1,2,2);                                    % ---- unemployment ----
fill([d; flipud(d)], [uA; flipud(uC)], [0.98 0.92 0.88], 'EdgeColor','none'); hold on  % drop on Octave
hA2 = plot(d, uA, '-',  'LineWidth', 3, 'Color', [0.15 0.15 0.15]);
hC2 = plot(d, uC, '--', 'LineWidth', 3, 'Color', [0.85 0.45 0.10]);
datetick('x','yyyy-mm','keeplimits'); xlim([d(1) d(end)]);
grid on; set(gca,'FontSize',12,'LineWidth',1.2);
ylabel('unemployment rate, %'); title('Unemployment');
legend([hA2 hC2], {'actual','tighter BoC (+100bp)'}, 'Location','northwest','FontSize',11);
gU = uC(end)-uA(end);
text(d(end), uC(end), sprintf('  %+.2f pp', gU), 'FontSize',12, ...
     'Color',[0.85 0.45 0.10], 'VerticalAlignment','middle');

lightTheme(gcf);
print(gcf, fullfile(outdir,'tighter_boc_price_unemp.png'), '-dpng', '-r150');

%% =====================================================================
%  Counterfactual 2 -- "what if the 2022 commodity spike hadn't happened?"
%  Same two-panel figure (price level + unemployment) as counterfactual 1.
%  =====================================================================
% oil enters in log-differences, so "hold the LEVEL at its Dec-2021 value"
% is a path of ZEROS in transformed units (no monthly growth), not r0.
C2 = wzCounterfactual(Rb.model, Rb.tr.Y, Rb.tr.dates, [2021 12], ...
                      ixo, zeros(Hcf,1), Hcf, [], ixo);   % confine to the oil shock
d2 = C2.dates;

pA2 = 100*cumsum(C2.actual(:,ixc));  pC2 = 100*cumsum(C2.cf(:,ixc));  % price: dlog -> cumulate
uA2 = C2.actual(:,ixu);              uC2 = C2.cf(:,ixu);              % unemp: already a level

fprintf('no 2022 commodity spike -> price level %+.2f pp vs actual after %d months\n', ...
        pC2(end)-pA2(end), Hcf);
fprintf('unemployment at end: actual %.2f%%, counterfactual %.2f%%\n', ...
        C2.actual(end,ixu), C2.cf(end,ixu));

figure('Color','w','Position',[100 100 1200 500]);

subplot(1,2,1);                                    % ---- price level ----
fill([d2; flipud(d2)], [pA2; flipud(pC2)], [0.90 0.95 0.90], 'EdgeColor','none'); hold on  % drop on Octave
hA = plot(d2, pA2, '-',  'LineWidth', 3, 'Color', [0.15 0.15 0.15]);
hC = plot(d2, pC2, '--', 'LineWidth', 3, 'Color', [0.15 0.55 0.25]);
datetick('x','yyyy-mm','keeplimits'); xlim([d2(1) d2(end)]);
grid on; set(gca,'FontSize',12,'LineWidth',1.2);
ylabel('cumulative price level, %'); title('Price level');
legend([hA hC], {'actual','no 2022 oil spike'}, 'Location','northwest','FontSize',11);
gP2 = pC2(end)-pA2(end);
text(d2(end), pC2(end), sprintf('  %+.2f pp', gP2), 'FontSize',12, ...
     'Color',[0.15 0.55 0.25], 'VerticalAlignment','middle');

subplot(1,2,2);                                    % ---- unemployment ----
fill([d2; flipud(d2)], [uA2; flipud(uC2)], [0.98 0.92 0.88], 'EdgeColor','none'); hold on  % drop on Octave
hA2 = plot(d2, uA2, '-',  'LineWidth', 3, 'Color', [0.15 0.15 0.15]);
hC2 = plot(d2, uC2, '--', 'LineWidth', 3, 'Color', [0.85 0.45 0.10]);
datetick('x','yyyy-mm','keeplimits'); xlim([d2(1) d2(end)]);
grid on; set(gca,'FontSize',12,'LineWidth',1.2);
ylabel('unemployment rate, %'); title('Unemployment');
legend([hA2 hC2], {'actual','no 2022 oil spike'}, 'Location','northwest','FontSize',11);
gU2 = uC2(end)-uA2(end);
text(d2(end), uC2(end), sprintf('  %+.2f pp', gU2), 'FontSize',12, ...
     'Color',[0.85 0.45 0.10], 'VerticalAlignment','middle');

lightTheme(gcf);
print(gcf, fullfile(outdir,'no_oil_spike_price_unemp.png'), '-dpng', '-r150');

fprintf('run_counterfactuals: done. Figures in ./figures/counterfactual.\n');
