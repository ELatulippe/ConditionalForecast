function run_all()
%RUN_ALL  Run every exercise end to end, then assemble the PDF report.
%
%   >> run_all
%
% Executes each run_*.m script in order, catching errors so one failed exercise
% does not stop the rest, times each, prints a summary, and finally calls
% build_report to collect the figures and the FEVD table into figures/report.pdf.
%
% OFFLINE: with the shipped panel caches present (cache/panel_k9.mat and
% cache/panel_k13.mat), no exercise contacts FRED -- main_risk_scenarios reads
% the matching cache and returns it. If a cache is missing, the affected steps
% will try to fetch from FRED (needs internet); build it once with freeze_panel
% to guarantee a fully offline, reproducible run.

    root = fileparts(mfilename('fullpath'));
    addpath(root); addpath(fullfile(root,'src'));
    assert(exist('main_risk_scenarios','file')==2, ...
        'main_risk_scenarios.m not found under ./src -- keep this script at the project root.');

    % Deterministic seed for any stochastic step (fan-chart / IRF draws).
    try, rng(20260727);
    catch, rand('state',20260727); randn('state',20260727);
    end

    % Cache pre-check -> tell the user whether this run is offline.
    ck9  = exist(fullfile(root,'cache','panel_k9.mat'),'file')  == 2;
    ck13 = exist(fullfile(root,'cache','panel_k13.mat'),'file') == 2;
    if ck9 && ck13
        fprintf('run_all: panel caches present -- running OFFLINE (no FRED calls).\n\n');
    else
        fprintf(2, ['run_all: one or both panel caches are missing, so some steps\n' ...
                    'will fetch from FRED (needs internet). Run freeze_panel once to\n' ...
                    'make the whole pipeline offline and reproducible.\n\n']);
    end

    steps = { ...
        'run_main', ...
        'run_single_scenarios', ...
        'run_joint_scenarios', ...
        'run_temporary_shocks', ...
        'run_soft_conditions', ...
        'run_counterfactuals', ...
        'run_impulse_responses', ...
        'run_variance_decomp' };

    status = cell(numel(steps),3);   % name | OK/FAIL | seconds or message
    for i = 1:numel(steps)
        name = steps{i};
        fprintf('\n========== %s ==========\n', name);
        t0 = tic;
        try
            run(fullfile(root,[name '.m']));
            status(i,:) = {name, 'OK', sprintf('%.1fs', toc(t0))};
        catch err
            status(i,:) = {name, 'FAIL', err.message};
            fprintf(2, '  %s FAILED: %s\n', name, err.message);
        end
        close all;   % keep figure handles from piling up across steps
    end

    % ---- summary ----
    fprintf('\n================= run_all summary =================\n');
    nfail = 0;
    for i = 1:numel(steps)
        if strcmp(status{i,2},'FAIL'), nfail = nfail + 1; end
        fprintf('  %-24s %-4s  %s\n', status{i,1}, status{i,2}, status{i,3});
    end
    fprintf('  %d of %d steps OK\n', numel(steps)-nfail, numel(steps));

    % ---- assemble the report from whatever figures were produced ----
    fprintf('\n========== build_report ==========\n');
    try
        build_report(fullfile(root,'figures'));
    catch err
        fprintf(2, '  build_report FAILED: %s\n', err.message);
    end

    fprintf('\nrun_all: done.\n');
end
