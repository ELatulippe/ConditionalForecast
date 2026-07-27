function scen = buildScenarios(tr, H, spec)
% BUILDSCENARIOS  Construct the constraint paths fed to wzConditional.
%
%   scen = buildScenarios(tr, H)          % the paper's four scenarios
%   scen = buildScenarios(tr, H, spec)    % your own
%
% Each path is expressed in the VAR's TRANSFORMED units, which is what the
% conditioning machinery consumes:
%   - dlog variables (oil, US IP, CPI, housing, FX): monthly log-differences
%   - level variables (unemployment, the policy rate): the level itself
% Use NaN for any horizon you want left FREE. A path that is NaN after month
% 12 constrains only the first year and lets the model determine the rest --
% wzConditional and wzBands both skip NaN entries when building R and r.
%
% SPEC
%   spec is a CELL ARRAY of structs, one per scenario. A cell array rather
%   than a struct array because the fields differ by type -- a 'ramp' needs
%   .to and .months, a 'policy' needs .step and .peak -- and MATLAB refuses
%   to build a struct array whose elements have different fields. A struct
%   array is still accepted when your scenarios happen to share a field set.
%
%   Each entry has fields:
%
%     .var     variable name ('oil','usip','cpi','house','unemp','rate','fx')
%              or a column index. Omit it and give .parts instead to
%              constrain SEVERAL variables at once -- see JOINT SCENARIOS.
%     .type    one of:
%                'ramp'   ramp the LEVEL from its last observed value to
%                         .to over .months, then hold. Fields: .to, .months
%                'hold'   hold the LEVEL at .at (default: last observed)
%                'replay' repeat the variable's own historical growth from
%                         .from = [year month] for .months months
%                'policy' the step-hold-ease path. Fields: .step, .peak,
%                         .hold, .floor
%                'level'  supply .path yourself, in LEVELS (converted to
%                         log-differences for dlog variables)
%                'pulse'  ramp the LEVEL to .to over .months, hold it there
%                         for .hold months, then leave the variable FREE.
%                         Fields: .to, .months, .hold (default 0)
%                'raw'    supply .path yourself, already in transformed units
%     .name    label for legends and filenames
%     .color   1x3 RGB (optional)
%     .desc    one-line description (optional)
%     .window  optional [h1 h2]; the constraint applies only to horizons
%              h1..h2 and every other horizon is set NaN, i.e. left free
%     .sd      optional scalar or H-vector: the standard deviation of a SOFT
%              condition, in the units of the path. Omit or 0 for a hard
%              condition. A soft condition says "around this path" rather
%              than "exactly this path", and the constrained variable then
%              gets a band of its own instead of a degenerate line. For the
%              policy rate or unemployment the units are percentage points,
%              so .sd = 0.5 reads as plus or minus about one point at 90%.
%
% Any field you omit falls back to the paper's value for that type.
%
% JOINT SCENARIOS
%   Give .parts, a cell array of sub-specs, to constrain several variables
%   within one scenario. Each sub-spec has its own .var, .type and shape
%   fields; .name, .color and .desc come from the outer struct.
%
%     struct('name','Stagflation', 'parts',{{ ...
%        struct('var','oil',  'type','ramp','to',110,'months',6), ...
%        struct('var','unemp','type','hold','at',8.0) }})
%
%   The conditioning side has always accepted several variables -- R and r
%   simply gain rows -- so nothing downstream changes. Watch ||eps||: joint
%   conditions are far more demanding than either leg alone, and a scenario
%   that needs many times a typical draw is one the model regards as close to
%   impossible.
%
% PARTIAL-HORIZON CONDITIONS
%   NaN leaves a horizon free, and there are two ways to get there. .window
%   blanks everything outside [h1 h2]. Type 'pulse' blanks everything after
%   the plateau, which is usually what you want for a temporary shock:
%
%     struct('var','oil','type','pulse','to',110,'months',3,'hold',3, ...
%            'name','Oil spike, then whatever follows')
%
%   For a log-difference variable, constraining months 1..n pins the LEVEL at
%   month n and lets it evolve freely from there -- so the model supplies the
%   persistence instead of you imposing a plateau. Note the asymmetry: a
%   .window that does not start at h=1 constrains growth rates without having
%   pinned the level they build on, which is rarely what is meant.
%
% EXAMPLES
%   spec = { ...
%     % Oil to $95 over 3 months instead of $120 over 7
%     struct('var','oil','type','ramp','to',95,'months',3, ...
%            'name','Oil to $95'), ...
%     % Unemployment held at 4.5% for the first year only, free thereafter
%     struct('var','unemp','type','hold','at',4.5,'window',[1 12], ...
%            'name','Tight labour market, 1yr'), ...
%     % A shallower tightening: +25bp/month to 5.5%, hold a year, ease to 3%
%     struct('var','rate','type','policy','step',0.25,'peak',5.5, ...
%            'hold',12,'floor',3,'name','Gradual tightening'), ...
%     % Replay the 2020 collapse in U.S. activity instead of 2008-09
%     struct('var','usip','type','replay','from',[2020 1], ...
%            'name','U.S. pandemic-style collapse'), ...
%     % An arbitrary hand-built oil path, in levels
%     struct('var','oil','type','level', ...
%            'path',[linspace(80,140,12) 140*ones(1,36)], ...
%            'name','Oil spike then plateau') };
%   R = main_risk_scenarios('fred', setfield(opts,'scenarios',spec));
%
% OUTPUT: 1xM struct array with fields .name, .cvar, .path (H x 1), .color,
% .desc, .type.

    idx = tr.idx;  L = tr.levels;  dts = tr.rawdates;

    if nargin < 3 || isempty(spec)
        spec = defaultSpec();
    end
    if isstruct(spec)                      % accept a struct array too
        spec = num2cell(spec);
    end
    if ~iscell(spec)
        error('buildScenarios:spec', ...
            'spec must be a cell array of structs (or a struct array).');
    end

    palette = [0.55 0.20 0.65; 0.85 0.33 0.10; 0.00 0.45 0.74; ...
               0.93 0.69 0.13; 0.30 0.60 0.30; 0.60 0.30 0.30];

    scen = struct('name',{},'cvar',{},'path',{},'sd',{},'color',{},'desc',{},'type',{});
    for m = 1:numel(spec)
        sp = spec{m};

        if isfield(sp,'parts') && ~isempty(sp.parts)
            parts = sp.parts;
            if ~iscell(parts), parts = num2cell(parts); end
            cv = zeros(1, numel(parts));
            pth = nan(H, numel(parts));
            sdm = nan(H, numel(parts));
            tys = cell(1, numel(parts));
            for q = 1:numel(parts)
                [cv(q), pth(:,q), tys{q}] = onePath(parts{q}, tr, H, m);
                sdm(:,q) = sdVector(parts{q}, H);
            end
            if numel(unique(cv)) < numel(cv)
                error('buildScenarios:dupvar', ...
                    ['Scenario %d constrains the same variable twice. Merge ' ...
                     'those parts into one path.'], m);
            end
            typ = strjoin(tys, '+');
        else
            [cv, pth, typ] = onePath(sp, tr, H, m);
            sdm = sdVector(sp, H);
        end

        n = numel(scen) + 1;
        scen(n).name  = getf(sp,'name',  sprintf('Scenario %d', m));
        scen(n).cvar  = cv;
        scen(n).path  = pth;
        scen(n).sd    = sdm;
        scen(n).color = getf(sp,'color', palette(1+mod(m-1,size(palette,1)),:));
        scen(n).desc  = getf(sp,'desc',  '');
        scen(n).type  = typ;
    end
end

% ======================================================================
function [v, path, typ] = onePath(sp, tr, H, m)
% Build one variable's constraint path from a (sub-)spec.
    idx = tr.idx;  L = tr.levels;  dts = tr.rawdates;
    if ~isfield(sp,'var') || ~isfield(sp,'type')
        error('buildScenarios:spec', ...
            'Scenario %d needs at least .var and .type.', m);
    end
    v = resolveVar(sp.var, idx);
    isLevelVar = strcmp(tr.tcode{v}, 'level');
    typ = lower(sp.type);

    switch typ
        case 'ramp'
            to     = getf(sp,'to',    120);
            months = getf(sp,'months',7);
            lev = rampHold(L(end,v), to, months, H);
            path = levelToUnits(lev, L(end,v), isLevelVar);

        case 'pulse'
            to     = getf(sp,'to',    120);
            months = getf(sp,'months',6);
            holdN  = getf(sp,'hold',  0);
            nCon   = min(months + holdN, H);
            lev    = rampHold(L(end,v), to, months, nCon);
            p      = levelToUnits(lev, L(end,v), isLevelVar);
            path   = padTo(p(1:nCon), H);       % everything after is FREE

        case 'hold'
            at   = getf(sp,'at', L(end, v));
            path = levelToUnits(repmat(at,H,1), L(end,v), isLevelVar);

        case 'replay'
            from   = getf(sp,'from',  [2008 1]);
            months = getf(sp,'months', H);
            path   = replayGrowth(L(:,v), dts, from, months, H, isLevelVar);

        case 'policy'
            if ~isLevelVar
                error('buildScenarios:policyOnGrowth', ...
                    ['type ''policy'' builds a level path, but %s enters ' ...
                     'the VAR in log-differences.'], tr.names{v});
            end
            path = policyPath(L(end,v), H, ...
                              getf(sp,'step',0.5), getf(sp,'peak',6), ...
                              getf(sp,'hold',8),   getf(sp,'floor',4));

        case 'level'
            lev  = padTo(sp.path(:), H);
            path = levelToUnits(lev, L(end,v), isLevelVar);

        case 'raw'
            path = padTo(sp.path(:), H);

        otherwise
            error('buildScenarios:type', 'Unknown scenario type "%s".', sp.type);
    end

    if isfield(sp,'window') && ~isempty(sp.window)
        w = sp.window;  keep = false(H,1);
        keep(max(1,w(1)) : min(H,w(2))) = true;
        path(~keep) = NaN;
    end
    path = path(:);
end

function sd = sdVector(sp, H)
% Per-horizon standard deviation of a soft condition; NaN means hard.
    if ~isfield(sp,'sd') || isempty(sp.sd)
        sd = nan(H,1);  return;
    end
    v = sp.sd(:);
    if numel(v) == 1, v = repmat(v, H, 1); end
    if numel(v) < H,  v = [v; nan(H-numel(v),1)]; end
    sd = v(1:H);
end

function lev = rampHold(x0, to, months, n)
    lev = zeros(n,1);
    for h = 1:n
        if h <= months, lev(h) = x0 + (to - x0)*(h/months);
        else,           lev(h) = to;
        end
    end
end

% ======================================================================
function spec = defaultSpec()
% The four scenarios of Section 4, as a spec.
    spec = { ...
      struct('var','oil','type','ramp','to',120,'months',7, ...
             'name','High oil prices','color',[0.55 0.20 0.65], ...
             'desc','Oil rises to $120 over 7 months then holds.'), ...
      struct('var','usip','type','replay','from',[2008 1], ...
             'name','U.S. recession','color',[0.85 0.33 0.10], ...
             'desc','US industrial production replays its 2008-09 growth path.'), ...
      struct('var','unemp','type','hold', ...
             'name','Low unemployment','color',[0.00 0.45 0.74], ...
             'desc','Unemployment held at its final observed value.'), ...
      struct('var','rate','type','policy','step',0.5,'peak',6, ...
             'hold',8,'floor',4, ...
             'name','Tight monetary policy','color',[0.93 0.69 0.13], ...
             'desc','Bank rate raised to 6%, held 8 months, eased to 4%.') };
end

% ======================================================================
function v = resolveVar(spec, idx)
    if ischar(spec)
        if ~isfield(idx, spec)
            error('buildScenarios:var', ...
                'Unknown variable "%s". Use one of: %s.', spec, ...
                strjoin(fieldnames(idx)', ', '));
        end
        v = idx.(spec);
    else
        v = spec;
    end
end

function y = getf(s, f, dflt)
    if isfield(s, f) && ~isempty(s.(f)), y = s.(f); else, y = dflt; end
end

function x = padTo(x, H)
% Shorter path -> the remaining horizons are left FREE, not extrapolated.
    if numel(x) < H
        x = [x; nan(H-numel(x), 1)];
    else
        x = x(1:H);
    end
end

function path = levelToUnits(lev, x0, isLevelVar)
% Convert a LEVEL path into the units the VAR uses for that variable. A NaN
% level yields a NaN growth rate on both sides of it, which is what "free"
% means -- the level at that point is not pinned, so neither is the step into
% or out of it.
    if isLevelVar
        path = lev(:);
    else
        path = diff(log([x0; lev(:)]));
    end
end

% ======================================================================
function path = replayGrowth(x, dts, from, months, H, isLevelVar)
% Repeat the variable's own historical behaviour starting just after `from`.
    j = find(dts == datenum(from(1), from(2), 1), 1);
    if isLevelVar
        src = x;                       % levels replay directly
        off = 0;
    else
        src = diff(log(x));            % element m is dated dts(m+1)
        off = 0;                       % so rows j..j+n-1 are the months AFTER `from`
    end
    n = min(months, H);
    if isempty(j) || j + off + n - 1 > numel(src)
        warning('buildScenarios:replay', ...
            ['History from %04d-%02d is too short to replay %d months; ' ...
             'using the last %d available instead.'], from(1), from(2), n, n);
        seg = src(end-n+1:end);
    else
        seg = src(j+off : j+off+n-1);
    end
    path = padTo(seg(:), H);           % horizons beyond `months` left free
end

% ======================================================================
function path = policyPath(r0, H, step, peak, holdN, floorR)
% Step up to a peak, pause, then ease to a floor. The month the peak is
% reached counts as the first month of the pause, so the plateau lasts
% exactly holdN months.
    path = zeros(H,1);  r = r0;  phase = 'up';  held = 0;
    for h = 1:H
        switch phase
            case 'up'
                r = min(r + step, peak);
                if r >= peak - 1e-12, phase = 'hold'; held = 0; end
            case 'hold'
                % rate unchanged during the pause
            case 'down'
                r = max(r - step, floorR);
        end
        path(h) = r;
        if strcmp(phase,'hold')
            held = held + 1;
            if held >= holdN, phase = 'down'; end
        end
    end
end
