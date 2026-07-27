function chk = checkTable1(shares, idx, verbose)
% CHECKTABLE1  Compare the computed FEV decomposition with the paper's Table 1.
%
%   chk = checkTable1(shares, idx)
%   chk = checkTable1(shares, idx, false)   % silent
%
% Table 1 depends ONLY on the estimated VAR -- not on any scenario -- so it
% is the cleanest test of whether the data layer is right. If it does not
% land, no conditional forecast will, and there is no point debugging the
% scenario figures yet.
%
% INPUTS
%   shares   K x K x Hmax from varianceDecomp; shares(k,l,h) is the share of
%            the h-step forecast-error variance of variable k due to shock l
%   idx      struct of column indices from loadCanadaData
%   verbose  print the side-by-side table (default true)
%
% OUTPUT struct 'chk'
%   .horizons  [3 12 24 48]
%   .computed  5 x 4 x 4  (variable x shock x horizon), in percent
%   .paper     5 x 4 x 4  the paper's point estimates
%   .inCI      5 x 4 x 4  logical, computed value inside the paper's 95%
%              bootstrap interval
%   .fracInCI  share of the 80 cells inside the interval
%   .gates     struct of the three headline checks
%   .pass      logical, all three gates satisfied
%
% Reference values are transcribed from Table 1, Appendix A. Note that the
% published interval for the h=3 response of unemployment to the U.S.
% activity shock ([0 3]) excludes its own point estimate (6); that cell's
% interval is widened here so the typo does not register as a failure.

    if nargin < 3 || isempty(verbose), verbose = true; end

    H    = [3 12 24 48];
    Hmax = size(shares, 3);
    avail = H <= Hmax;
    if ~any(avail)
        error('checkTable1:horizon', ...
            ['The decomposition only reaches h=%d; Table 1 starts at h=3. ' ...
             'Raise opts.H.'], Hmax);
    end
    if ~all(avail)
        fprintf(2, ['Horizons %s exceed the forecast horizon (H=%d) and are ' ...
                    'dropped from the\nTable 1 comparison.\n'], ...
                mat2str(H(~avail)), Hmax);
    end
    vrow = [idx.cpi idx.house idx.unemp idx.rate idx.fx];
    scol = [idx.oil idx.usip idx.unemp idx.rate];
    vname = {'Inflation','Housing','Unemployment','Bank rate','Exch. rate'};
    sname = {'oil','US activity','unemployment','monetary policy'};

    [paper, lo, hi] = table1Reference();
    H = H(avail);  paper = paper(:,:,avail);  lo = lo(:,:,avail);  hi = hi(:,:,avail);

    computed = zeros(5,4,numel(H));
    for ih = 1:numel(H)
        computed(:,:,ih) = 100 * shares(vrow, scol, H(ih));
    end
    inCI = (computed >= lo - 1e-12) & (computed <= hi + 1e-12);

    % ---- the three headline gates ------------------------------------
    oilInfl  = squeeze(computed(1,1,:));    % oil -> inflation
    oilUnemp = squeeze(computed(3,1,:));    % oil -> unemployment
    mpRate   = squeeze(computed(4,4,:));    % MP  -> bank rate
    gates = struct( ...
        'oilOnInflation',   struct('value',oilInfl,  'target',[44 60],'ok',all(oilInfl  >= 24 & oilInfl  <= 70)), ...
        'oilOnUnemployment',struct('value',oilUnemp, 'target',[16 22],'ok',all(oilUnemp >=  5 & oilUnemp <= 41)), ...
        'mpOnBankRate',     struct('value',mpRate,   'target',[62 96],'ok',all(mpRate   >= 31 & mpRate   <= 98)));
    pass = gates.oilOnInflation.ok && gates.oilOnUnemployment.ok && gates.mpOnBankRate.ok;

    chk = struct('horizons',H,'computed',computed,'paper',paper, ...
                 'inCI',inCI,'fracInCI',mean(inCI(:)),'gates',gates,'pass',pass);

    if ~verbose, return; end

    fprintf('\n================ Table 1 check (computed vs paper) ================\n');
    fprintf('%% of forecast-error variance; "*" = outside the paper''s 95%% interval\n');
    for s = 1:4
        fprintf('\nShock to %s\n', sname{s});
        fprintf('  %-14s', 'variable');
        for ih = 1:numel(H), fprintf('%16s', sprintf('h=%d', H(ih))); end
        fprintf('\n');
        for v = 1:5
            fprintf('  %-14s', vname{v});
            for ih = 1:numel(H)
                flag = ' '; if ~inCI(v,s,ih), flag = '*'; end
                fprintf('%12s%c   ', sprintf('%.0f (%.0f)', computed(v,s,ih), paper(v,s,ih)), flag);
            end
            fprintf('\n');
        end
    end
    fprintf('\n  cells inside the paper''s 95%% intervals: %d / %d (%.0f%%)\n', ...
            sum(inCI(:)), numel(inCI), 100*chk.fracInCI);

    fprintf('\n  Headline gates\n');
    printGate('oil -> inflation      ', oilInfl,  '44-60%');
    printGate('oil -> unemployment   ', oilUnemp, '16-22%');
    printGate('MP  -> bank rate      ', mpRate,   '62-96%');
    if pass
        fprintf('\n  GATES PASSED -- the data layer looks right; scenario figures are worth reading.\n');
    else
        % Name the specific cells that fail rather than reciting a fixed
        % checklist, which goes stale as each cause is ruled out.
        fprintf(2, '\n  GATES FAILED at:\n');
        gnames = {'oil -> inflation','oil -> unemployment','MP  -> bank rate'};
        gvals  = {oilInfl, oilUnemp, mpRate};
        glo    = [24  5 31];  ghi = [70 41 99];
        for g = 1:3
            bad = find(gvals{g} < glo(g) | gvals{g} > ghi(g));
            if isempty(bad), continue; end
            fprintf(2, '    %-22s', gnames{g});
            for b = bad(:).'
                fprintf(2, ' h=%d (%.0f vs %.0f)', H(b), gvals{g}(b), paper(gaterow(g), gatecol(g), b));
            end
            fprintf(2, '\n');
        end
        fprintf(2, ['  Usual causes, roughly in order: the house-price series, the terminal\n' ...
                    '  month of the panel, seasonal adjustment of CPI, the lag order, and the\n' ...
                    '  2020 episode. Once those are ruled out, a gate that fails only at short\n' ...
                    '  horizons points at Sigma_u; one that fails only at long horizons points\n' ...
                    '  at the dynamics.\n']);
    end
    fprintf('==================================================================\n\n');
end

% ----------------------------------------------------------------------
function r = gaterow(g)
    rows = [1 3 4];  r = rows(g);      % inflation, unemployment, bank rate
end

function c = gatecol(g)
    cols = [1 1 4];  c = cols(g);      % oil, oil, monetary policy
end

% ----------------------------------------------------------------------
function printGate(label, v, target)
    fprintf('    %s paper %-8s computed %s\n', label, target, ...
            sprintf('%.0f ', v));
end

% ----------------------------------------------------------------------
function [pt, lo, hi] = table1Reference()
% Table 1 of Appendix A, laid out as (variable x shock x horizon).
% Variables: inflation, housing, unemployment, bank rate, exchange rate.
% Shocks:    oil, US activity, unemployment, monetary policy.
% Horizons:  3, 12, 24, 48 months.

    % Panel A: oil                     h=3          h=12         h=24         h=48
    A  = [44 60 49 45;   5  3  5  5;  19 22 17 16;   2  3  3  3;  20 14 14 14];
    Al = [33 38 27 24;   1  1  2  2;   5  5  5  5;   0  0  0  1;  10  7  7  7];
    Ah = [55 70 62 59;  12 14 18 20;  34 41 36 34;   8 17 20 18;  31 27 28 27];

    % Panel B: U.S. economic activity
    B  = [ 1  1  1  4;   0  2  8 11;   6  7  6  9;   1  7 10  9;   1  2  3  3];
    Bl = [ 0  0  0  1;   0  0  0  1;   0  0  1  2;   0  0  1  1;   0  0  0  1];
    Bh = [ 5  7 10 14;   4 12 24 27;   6 25 22 23;   6 21 29 30;   5 12 16 15];
    %                                  ^ published as [0 3]; widened to cover
    %                                    its own point estimate of 6 (typo).

    % Panel C: Canadian unemployment
    C  = [ 0  1  4  4;   0  2  2  3;  74 66 52 46;   0  1  8 16;   1  1  1  2];
    Cl = [ 0  0  0  1;   0  0  0  1;  46 34 27 23;   0  0  0  1;   0  0  0  1];
    Ch = [ 1  8 16 17;   2 10 15 18;  90 80 65 59;   4 10 25 35;   4  8 10 12];

    % Panel D: Canadian monetary policy
    D  = [ 1  1  2  3;   0  0  1  1;   1  1  5  7;  96 87 71 62;   0  0  1  1];
    Dl = [ 0  0  0  1;   0  0  0  0;   0  0  1  1;  85 63 41 31;   0  0  0  0];
    Dh = [ 3  6 13 15;   1  6 15 17;   5  8 17 22;  98 91 82 77;   2  6 10 12];

    pt = zeros(5,4,4);  lo = pt;  hi = pt;
    P = {A,B,C,D};  Lo = {Al,Bl,Cl,Dl};  Hi = {Ah,Bh,Ch,Dh};
    for s = 1:4
        for ih = 1:4
            pt(:,s,ih) = P{s}(:,ih);
            lo(:,s,ih) = Lo{s}(:,ih);
            hi(:,s,ih) = Hi{s}(:,ih);
        end
    end
end
