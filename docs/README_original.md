# Risk Scenarios and Macroeconomic Forecasts — MATLAB replication

MATLAB implementation of the conditional-forecast methodology in
Moran, Stevanovic & Surprenant (2025, BoC SWP 2025-28), following the
equation-by-equation derivation note. It estimates a monthly VAR for the
Canadian economy, produces a baseline forecast, and computes four risk
scenarios (high oil prices, U.S. recession, low unemployment, tight monetary
policy) via the Waggoner–Zha (WZ) method, with a Baumeister–Kilian (BK)
robustness check and a forecast-error variance decomposition.

Runs on **base MATLAB** (R2014b+) — no toolboxes required (OLS, Cholesky, and
the min-norm solve are hand-rolled). Also runs in **GNU Octave** (tested on 8.4).

## Quick start

You need one file the package cannot download for you: a **monthly** Canadian
house-price index (see "Housing" below).

```matlab
addpath(pwd)
opts = struct('houseFile', '1810020501-eng.csv');   % StatCan NHPI, as downloaded
R = main_risk_scenarios('fred', opts);          % pull the rest from FRED

% offline, from CSVs you already downloaded:
% R = main_risk_scenarios('local', struct('csvDir','./csv', ...
%                                         'houseFile','1810020501-eng.csv'));
% R = main_risk_scenarios('synthetic');          % fake data, code test only
```

Figures are written to `./figures/`. `R` holds the model, forecasts, scenarios,
the variance decomposition, and `R.table1` (the Table 1 check).

Validate the numerical engines independently with `test_engines`.

## What the driver does, in order

The run is sequenced so that the things most likely to be wrong are caught
before anything downstream is computed:

1. **Load the panel and assert the terminal month.** Every scenario reads its
   starting level off the last row (`oil0`, `u0`, `r0`), so a silently trimmed
   sample shifts the whole exercise. `loadCanadaData` now prints a per-series
   first/last table and names whichever series binds; `main_risk_scenarios`
   errors out if the panel does not reach the requested end date. Override
   with `opts.requireEnd = false`.
2. **Estimate at the paper's lag order.** `p = 3` is fixed, not selected. AIC
   is still computed and reported as a diagnostic — if it disagrees, that is
   usually a sign that CPI is not seasonally adjusted and the criterion is
   buying lags to chase the seasonal pattern.
3. **Print the variance decomposition and check it against Table 1.**
   `checkTable1` compares all 80 cells against the published point estimates
   and 95% bootstrap intervals, and evaluates three headline gates:
   oil → inflation (44–60%), oil → unemployment (16–22%), and monetary policy
   → bank rate (62–96%). Table 1 depends only on the estimated VAR, so it is
   the cleanest test of the data layer.
4. **Scenarios and figures** — skipped by default when the gate fails, since
   there is nothing to learn from the figures until the VAR itself matches.
   Draw them anyway with `opts.gate = false`.

## Options

| Field | Default | Meaning |
|---|---|---|
| `houseFile` | *(required)* | monthly house-price CSV (long or StatCan wide) |
| `houseSeries` | `Total (house and land)` | which row to take from a wide export |
| `cpiFile` | `''` | seasonally adjusted CPI CSV, overriding the FRED default |
| `cpiSeries` | `All-items` | which row to take from a wide CPI export |
| `extra` | `{}` | add variables to the panel (see below) |
| `ids` | `struct()` | override the FRED series ID for any variable |
| `files` | `struct()` | CSV path for ANY series, keyed by variable name |
| `rows` | `struct()` | wide-export row label, keyed the same way |
| `rateFile` | `''` | monthly policy/short-rate CSV, overriding the FRED Bank Rate |
| `rateSeries` | `''` | which row to take from a wide rate export |
| `csvDir` | `pwd` | directory of pre-downloaded FRED CSVs (`source = 'local'`) |
| `cacheFile` | `''` | `.mat` file to save/reload the assembled panel |
| `allowHouseSpline` | `false` | fall back to the quarterly BIS proxy (testing only) |
| `p` | `3` | lag order |
| `H` | `48` | forecast horizon, months |
| `sample` | `[1992 1 2022 12]` | estimation window; give `[1992 1]` to run to the latest available month |
| `requireEnd` | `true` | error if the panel does not reach the terminal month |
| `covid` | `[]` | pandemic dummies: `'default'` = Mar–Aug 2020, or `[y1 m1; y2 m2]` |
| `seasonalDummies` | `false` | append 11 month dummies to the exogenous block |
| `prior` | `[]` | Minnesota prior, e.g. `struct('lambda',0.2)` |
| `gate` | `true` | skip figures when the Table 1 check fails |
| `figures` | `true` | draw figures at all |
| `paperFigs` | `true` | draw the paper-replication figures 1–4 |
| `diffFigs` | `true` | draw the scenario-minus-baseline figures 6 |
| `figStart` | `[2019 1]` | year and month where the figure x-axis begins |
| `figEnd` | `[]` | year and month where it ends (default: the full horizon) |
| `scenarios` | `[]` | scenario spec (default: the paper's four) — see `buildScenarios` |
| `bands` | `false` | compute uncertainty bands and draw a fan chart per scenario |
| `bandOpts` | `struct()` | passed to `wzBands` (`nDraws`, `parameter`, `shock`, ...) |

## Data map

| # | Variable | Source | Transform |
|---|---|---|---|
| 1 | WTI crude, USD/bbl | FRED `MCOILWTICO` | dlog |
| 2 | US industrial production | FRED `INDPRO` | dlog |
| 3 | Canada CPI, all items | FRED `CPALCY01CAM661N` (**NSA**) | dlog |
| 4 | Canada house-price index | **you supply** (`opts.houseFile`) | dlog |
| 5 | Canada unemployment rate | FRED `LRUNTTTTCAM156S` | level |
| 6 | Canada 3-month interbank rate | FRED `IR3TIB01CAM156N` | level |
| 7 | CAD per USD | FRED `EXCAUS` | dlog |

Variable **ordering matters** for the recursive (Cholesky) identification used
by BK and the variance decomposition: global/US variables first, sticky prices
before the activity/policy block, forward-looking FX last.

### Housing

FRED carries no monthly Canadian house-price index. The package therefore
**requires** you to supply one rather than silently substituting a proxy.

The paper's series is the monthly index in the Fortin-Gagnon et al. (2022)
database. Its 12-month growth is roughly 0% in 2019, peaks near 11.5% in late
2021, and sits around 4% at the end of 2022 — the profile in Figure 2. That
matches Statistics Canada's **New Housing Price Index** (table 18-10-0205); it
does *not* match the BIS residential property price index, whose 12-month
growth peaks near 26% over the same window.

`readSeriesCSV` also reads Bank of Canada multi-column exports, whose
observation block carries every sub-index side by side under a metadata
preamble. Name the column you want as the third argument, or as `opts.rows`:

```matlab
ob.files = struct('oil','BCPI_MONTHLY.csv');
ob.rows  = struct('oil','M.ENER');       % energy sub-index
```

With more than two value columns and no name given it errors with the list
rather than defaulting to the second column — silently returning the total
index in place of the energy one would be hard to notice.

`readSeriesCSV` accepts the StatCan **"Download as displayed"** export
directly — the wide layout with months running across a header row and one
series per row beneath it — as well as a plain long `date,value` file. The
leading metadata block, symbol legend, table corrections and footnotes are
skipped, data-quality flags appended to a value (`40.6E`) are stripped, and
English and French month names are both recognised, so `-eng` and `-fra`
downloads work interchangeably. By default the `Total (house and land)` row is
taken; select another with `opts.houseSeries` (an unrecognised label errors
out with the list of rows the file actually contains).

Sanity check on the resulting column: 12-month growth should be roughly −0.2%
in mid-2019, peak near 12% in August 2021, and sit at about 3.9% in December
2022, with the December 2016 base equal to 100.

Setting `allowHouseSpline = true` restores the old behaviour (quarterly BIS,
spline-interpolated to monthly) behind a loud warning. Do not use it to
replicate: monthly growth from a splined quarterly index is about **95%
predictable from its own three lags**, with an innovation standard deviation
roughly one-fifth of the series' own. Because housing sits fourth in the
Cholesky ordering, that collapsed innovation variance contaminates `Sigma_u`,
`D`, every `Theta_j`, the variance decomposition, and all four conditional
forecasts — not just the housing panel.

### Seasonality

`CPALCY01CAM661N` is not seasonally adjusted, while the paper's database is.
Because a 12-month growth rate nets out deterministic seasonality, the
*figures* look plausible while the estimated monthly dynamics are wrong — so
`seasonalityCheck` runs on the estimation-frequency data straight after the
transform, regressing each growth rate on eleven month dummies and reporting
the R-squared and joint F.

This matters for the decomposition, not just for tidiness. A deterministic
seasonal is, from the VAR's point of view, a large part of that variable's
innovation that is orthogonal to every other variable. It inflates that
variable's OWN share of its forecast-error variance and deflates every other
shock's — the signature to look for when one row of Table 1 is too
self-explained.

The fix is Statistics Canada table 18-10-0006-01, *Consumer Price Index,
monthly, seasonally adjusted*, which downloads in the same wide layout:

```matlab
opts.cpiFile = '1810000601-eng.csv';   % row defaults to 'All-items'
```

For reference, the NHPI has no meaningful deterministic seasonality
(R-squared about 0.03 on the month dummies), so it needs no adjustment.

## Display conventions

The VAR is estimated on monthly log-differences; these are display-only
transforms (see `dispSeries.m`):

- Prices, housing, U.S. industrial production, exchange rate: **12-month
  percent change**, matching the paper's axes. Annualised monthly growth
  (`×1200`) was used previously for U.S. IP and the exchange rate and cannot
  match — April 2020 U.S. industrial production alone annualises to about
  −160%, while the paper's panel runs over ±20.
- Unemployment, policy rate: levels.
- Housing and FX (Fig. 3): reconstructed index levels.

## Minnesota prior

At K = 9 and p = 12 the VAR has 123 regressors per equation on roughly 400
months — about 3.3 observations per parameter, which OLS cannot support.
`opts.prior` shrinks instead:

```matlab
ob.p     = 12;
ob.prior = struct('lambda', 0.2);     % smaller = tighter
```

Implemented the way Banbura, Giannone and Reichlin do it: as dummy
observations appended to the sample, so the posterior mean is OLS on the
augmented data and the estimator itself is unchanged. Each equation is shrunk
toward a univariate random walk (variables in levels) or white noise
(variables already differenced), with longer lags shrunk harder. `lambda` sets
overall tightness — as it grows the estimates return to OLS exactly, which
`test_engines` checks.

Exogenous regressors get all-zero dummy columns, i.e. a flat prior: pandemic
and seasonal dummies should not be shrunk toward anything.

**Uncertainty is then sampled, not bootstrapped.** With a conjugate prior the
posterior is normal-inverse-Wishart in closed form, so both `irfCholesky` and
`wzBands` draw from it via `varPosteriorDraw`. Resampling residuals around a
shrunk point estimate would ignore the prior's own contribution to posterior
uncertainty, understating it exactly where the prior is doing the work.
`IRF.bandSource` and `B.bandSource` record which was used; `forceBoot = true`
overrides, and the bootstrap then reapplies the prior to every draw so it
resamples the same estimator. Kilian's bias adjustment is skipped under a
prior — it corrects an OLS small-sample bias that shrinkage has displaced.

The prior also largely cures the explosive-draw problem: on a test panel the
rejection rate fell from over half to about 15%, so the bands stop being
heavily conditional on stationarity.

**Choosing lambda.** 0.2 is a reasonable default but arbitrary. Report a
sweep, or select it the way Banbura, Giannone and Reichlin do, by matching
in-sample fit to a small benchmark VAR:

```matlab
for lam = [0.05 0.1 0.2 0.5 1]
    o = ob;  o.prior = struct('lambda', lam);
    Rl = main_risk_scenarios('fred', o);
    I  = irfCholesky(Rl, 48, struct('scaleTo',1,'verbose',false));
    fprintf('lambda %-5.2g  max|eig| %.4f  cpi h=12 %+.3f  h=24 %+.3f\n', ...
            lam, max(abs(eig(Rl.model.A))), ...
            I.resp(Rl.data.idx.cpi, Rl.data.idx.rate, 12), ...
            I.resp(Rl.data.idx.cpi, Rl.data.idx.rate, 24));
end
```

## Impulse responses

```matlab
IRF = irfCholesky(R, 48, struct('scaleTo',1,'boot',500));
figureIRF(IRF, {'rate','oil'}, 'figures');
```

Two choices decide whether the plot is readable.

**Cumulation.** Log-difference variables respond in GROWTH, but a price puzzle
is a statement about the price LEVEL, so those responses are cumulated by
default and read as percent deviations of the level. Level variables read in
percentage points.

**Normalisation.** A raw response is per one standard deviation of the shock,
and if the shock is small the response looks small for reasons unrelated to
transmission — a smooth market rate can have a policy innovation of a few
basis points. `scaleTo = 1` rescales so the shocked variable moves by one unit
on impact, giving responses per 100bp for a rate. `IRF.scale` reports the
unscaled impact move; read it before judging any magnitude.

Bands use the same residual bootstrap as `wzBands`, with each draw
renormalised by its own impact scale so a 100bp band really is a 100bp band.

## Which norm the scenario minimises

`R s = r` has infinitely many solutions and the scenario picks one by
minimising a norm. Doan, Litterman and Sims (1984), whom the paper invokes at
this step, minimise the sum of squares of the **structural** shocks. That is
the "most likely combination of shocks" reading: under `u ~ N(0, Sigma_u)` the
most likely `u` satisfying the constraint minimises `u' Sigma_u^{-1} u`,
equivalently `eps'eps` with `u = D eps`.

Minimising the plain Euclidean norm of the **reduced-form** innovations
instead treats every equation as equally expensive to shock. In this system
the innovation standard deviations span two orders of magnitude — oil around
0.09 in logs against CPI around 0.003 — so the unweighted solution buys
cheap-looking reductions in `||u||` by loading shocks that are enormous in
standard-deviation terms onto the low-variance equations.

Both solutions satisfy `R s = r` exactly, so the constrained variable follows
its path either way and the error is invisible in the panel you constrained.
It surfaces everywhere else. On the tight-monetary-policy scenario the
unweighted solution needs `||eps|| = 147` against 22 for the correct one, and
the implied deviations from baseline are inflated by a factor of 17 for CPI
and 42 for housing.

`wzConditional` defaults to `'structural'`. Pass `'reduced'` as the seventh
argument to reproduce the old behaviour for comparison.

## The 2020 episode

In April–May 2020 the oil price, U.S. industrial production and the Canadian
unemployment rate all moved by many times their historical standard
deviations in the same months, and those observations dominate `Sigma_u`.
Because the recursive ordering puts oil first and U.S. activity second, the
common pandemic movement is attributed to those two shocks: their share of
the unemployment forecast-error variance is inflated, unemployment's own
share is pushed to the floor, and the oil-to-CPI pass-through is flattened
because a 57% monthly oil innovation is not matched proportionally in
consumer prices.

Two ways to find out whether this is what separates your estimates from
Table 1. Re-estimate before the pandemic:

```matlab
o = opts; o.sample = [1992 1 2019 12]; o.gate = false; o.figures = false;
R19 = main_risk_scenarios('fred', o);
```

or absorb the episode with dummies while keeping the full sample:

```matlab
o = opts; o.covid = 'default';    % March-August 2020
Rc = main_risk_scenarios('fred', o);
```

Both are diagnostics. The pre-pandemic window moves the scenario start date,
and the paper uses no dummies, so neither should be presented as a
reproduction — they exist to localise the discrepancy.

## Real-time forecasts instead of replication

Pass a two-element `sample` to run to the latest month the data support:

```matlab
ob.sample    = [1992 1];        % no end date
ob.cacheFile = 'panel_live.mat';
```

Three things change automatically. The terminal-date assertion is skipped,
since there is nothing to assert against; the loader reports which series caps
the panel rather than warning about misalignment; and the Table 1 gate is
switched off, because the paper's decomposition was estimated on its own
window and is not a valid pass/fail target for a different one. The comparison
is still printed as a diagnostic.

The loader prints the last observed month for every series before it trims,
followed by `Panel capped at ... by: ...`. Any FRED series that stops updating
will cap the panel, so read those lines first. Any of the seven can then be
supplied from a CSV in either supported layout:

```matlab
opts.files = struct('unemp','1410028701-eng.csv');   % StatCan LFS
opts.rows  = struct('unemp','Unemployment rate');
```

If the replacement is a StatCan table, `fetchStatCan` downloads the full-table
zip and writes the two-column CSV for you:

```matlab
fetchStatCan('10-10-0144-01')                                    % see what's inside
fetchStatCan('10100144', {'Rates','Bank rate'}, 'rate.csv')      % extract one series
opts.files = struct('rate','rate.csv');
```

Called without an output name it lists the distinct values of every label
column, so you can see what to filter on before committing. If the filters
still leave more than one series it says so and names the columns you have not
constrained yet, rather than silently taking the first match.

If the replacement also lives on FRED, overriding the ID is easier than
downloading anything:

```matlab
opts.ids = struct('unemp','LRUN64TTCAM156S');
```

`fetchFRED` requests `fredgraph.csv?id=<ID>` with no date bounds, so it always
pulls the full series — a short panel means the series itself has stopped, not
that the download was truncated. The OECD-sourced defaults for unemployment
(`LRUNTTTTCAM156S`) and the policy rate (`IRSTCB01CAM156N`) are the ones that
go stale; the WTI, U.S. industrial production and exchange-rate series are
maintained by other providers and keep updating.

Use a fresh `cacheFile`. The cache is read before the date range is applied,
so an existing `panel_sa.mat` would hand back the old 2022 panel; the driver
now warns when the cached panel falls short of what was requested.

Everything downstream follows the new terminal month on its own — the
scenarios read their starting levels off the last row, so the oil ramp starts
from the latest price and the policy path from the latest rate. Check the
`policy` scenario's `peak` after moving the sample: stepping to 6% from a
starting rate of 2.25% is a very different assumption from stepping there from
4.5%, and the `||eps||` figure will tell you how extreme the model finds it.
`figStart` is worth moving forward too, since the default 2019 start leaves a
long history on the axis.

## Horizon versus view

`opts.H` (default 48) is the forecast horizon and drives everything: scenario
path length, the rows of `R`, the MA coefficients, the bands. `opts.figEnd`
crops the plot only.

They are not interchangeable. Shortening `H` is a different exercise —
`R` loses rows, so the scenario imposes fewer restrictions, `||eps||` falls,
and the scenario looks more plausible simply because less is being asked of
the model. Cropping with `figEnd` leaves every estimate identical.

```matlab
ob.H       = 48;          % estimate over four years, as the paper does
ob.figEnd  = [2028 6];    % but plot only two
```

Note also that shortening `H` truncates the scenario shapes rather than
rescaling them: the default policy path needs about 20 months to step up,
pause and ease back, so at `H = 12` you would never reach the easing leg.

## Adding variables

`opts.extra` splices new variables into the panel. Each entry needs a name, a
source (`.id` for FRED or `.file` for a CSV, with `.row` for a wide export), a
transform, and a position in the recursive ordering:

```matlab
ob.extra = { ...
  struct('name','spread','file','termspread.csv','tcode','level','after','rate'), ...
  struct('name','tsx',   'file','tsx.csv',       'tcode','dlog', 'after','spread') };
```

`data.idx` gains a field per variable, so scenarios (`'var','spread'`),
the decomposition and the seasonality check all pick them up automatically.
The five hard-coded figure panels still show the paper's variables; read new
ones from `R.shares` and `R.Fscen`.

**Two things to weigh before adding anything.** Ordering is identification,
not bookkeeping: a variable at position k is assumed not to respond within the
month to anything after it. Financial variables that price continuously
generally belong late, near the exchange rate. And each variable costs `K*p+1`
more coefficients per equation *and* adds an equation — at K=7, p=3 you are
already estimating 22 coefficients per equation on ~410 months, so past nine
or so variables OLS gets thin and a Minnesota prior (Bańbura, Giannone and
Reichlin, which the paper already cites) is the honest next step. Adding
variables also dilutes every share in the decomposition, making Table 1
comparisons harder to read.

## Defining your own scenarios

`buildScenarios` is spec-driven. Pass a cell array of structs as
`opts.scenarios` — a cell array rather than a struct array because the fields
differ by type and MATLAB will not build a struct array whose elements have
different fields.

```matlab
spec = { ...
  struct('var','oil',  'type','ramp',  'to',95,'months',3, 'name','Oil to $95'), ...
  struct('var','unemp','type','hold',  'at',4.5,'window',[1 12], ...
                                       'name','Tight labour market, 1yr'), ...
  struct('var','rate', 'type','policy','step',0.25,'peak',5.5,'hold',12, ...
                                       'floor',3,'name','Gradual tightening'), ...
  struct('var','usip', 'type','replay','from',[2020 1], ...
                                       'name','U.S. pandemic-style collapse') };
opts.scenarios = spec;
```

Types: `ramp` (level to `.to` over `.months`, then hold), `pulse` (ramp, hold
`.hold` months, then leave the variable free), `hold` (level at `.at`),
`replay` (repeat the variable's own growth from `.from`), `policy` (step up by
`.step` to `.peak`, pause `.hold` months, ease to `.floor`), `level` (your own
path in levels) and `raw` (your own path already in the VAR's transformed
units).

**Joint scenarios.** Give `.parts` instead of `.var` to constrain several
variables at once:

```matlab
struct('name','Stagflation', 'parts',{{ ...
   struct('var','oil',  'type','ramp','to',110,'months',6), ...
   struct('var','unemp','type','hold','at',8.0) }})

struct('name','U.S. recession, no policy response', 'parts',{{ ...
   struct('var','usip','type','replay','from',[2008 1],'months',24), ...
   struct('var','rate','type','hold') }})
```

`R` and `r` simply gain rows, so nothing downstream changes. Joint conditions
are much more demanding than either leg alone — watch `||eps||`.

**Temporary shocks.** `pulse` constrains the ramp and plateau and then leaves
the variable free, so the model supplies the persistence instead of you
imposing a plateau to the horizon:

```matlab
struct('var','oil','type','pulse','to',110,'months',3,'hold',3, ...
       'name','Oil spike, 6 months')
```

For a log-difference variable, constraining months 1..n pins the LEVEL at
month n and lets it evolve from there. Note the asymmetry with `.window`: a
window that does not start at h=1 constrains growth rates without having
pinned the level they build on, which is rarely what is meant.

**Horizons.** `opts.H` sets the forecast length. To constrain only part of it,
use `.window` — `[1 12]` constrains the first year and leaves the rest NaN,
which `wzConditional` and `wzBands` skip when building `R` and `r`, so the
model determines those months. A `level` or `raw` path shorter than `H` is
padded with NaN rather than extrapolated, so it behaves the same way.

**Units.** Paths for `unemp` and `rate` are levels; every other variable
enters in log-differences, and `ramp`, `hold` and `level` convert for you.
Only `raw` expects transformed units.

**Sanity check.** Watch the `||eps||` figure the driver prints. It reports the
scenario in units of a typical draw, so anything much above 1x is a scenario
the model considers very unlikely — useful when inventing paths by hand.

## Soft conditions

A hard condition pins the constrained variable exactly, which is a strange
thing to do to something labelled a risk scenario. Waggoner and Zha also allow
conditioning on a distribution. Write the restriction as a noisy observation,
`r = R eps + e` with `e ~ N(0, Omega)`, and update the prior `eps ~ N(0, I)`:

    eps | r  ~  N( R'(RR' + Omega)^{-1} r ,  I - R'(RR' + Omega)^{-1} R )

`Omega = 0` recovers the hard case exactly, so both share one code path. Add
`.sd` to a scenario spec, in the units of the path:

```matlab
struct('var','rate','type','policy','peak',6,'sd',0.4, ...
       'name','Tightening, roughly')
```

For the policy rate and unemployment those units are percentage points. The
constrained variable then gets a band of its own rather than a degenerate
line, and `||eps||` falls, since the model no longer has to force an exact
path. `info.constraintResidual` is NaN under a soft condition — it is not
meant to bind.

## In-sample counterfactuals

Nothing restricts the algebra to the forecast horizon. Over a window that
actually happened, `y_actual = Fbase + M u_actual`, so a counterfactual asks
for a different `r` and the natural answer keeps every realised shock the
counterfactual does not require you to change:

    delta = pinv(R) (r_cf - r_actual),   y_cf = y_actual + M delta

The baseline drops out. This perturbs HISTORY rather than a forecast, so the
counterfactual inherits every realised shock outside the span of the
constraint — the 2022 oil shock and so on are all still there.

```matlab
% what if the Bank had held the rate at 0.25% through 2022-23?
C = wzCounterfactual(R.model, R.tr.Y, R.tr.dates, [2021 12], ...
                     R.tr.idx.rate, 0.25*ones(24,1), 24);
plot(C.dates, [C.actual(:,R.tr.idx.cpi), C.cf(:,R.tr.idx.cpi)])
```

Feeding history back as the target returns history exactly, which is the
cleanest test that the perturbation is minimal; `test_engines` checks it.

**Which shocks do the work.** Left unrestricted the perturbation spreads
across every structural shock: it finds the most likely way for the rate to
have taken a different path, not the most likely *policy* reason. Since the
rate is ordered late, the cheapest explanation for "the rate stayed low" is
often "the economy was weaker", so the solution loads on demand shocks and the
counterfactual returns a LOWER price level rather than a higher one. Read
`C.epsNorm` to see which shocks were used. To confine the perturbation to
policy shocks, pass the shock index as a ninth argument:

```matlab
C = wzCounterfactual(R.model, R.tr.Y, R.tr.dates, [2021 12], ...
                     R.tr.idx.rate, 0.25*ones(24,1), 24, [], R.tr.idx.rate);
```

With H constraints and H free policy shocks the system is square, so the path
is still hit exactly; `||delta||` rises, because the cheapest shocks are no
longer available. This is the Waggoner-Zha versus Baumeister-Kilian
distinction the paper raises in Section 6, applied to a counterfactual.

**The caveat belongs in any write-up using this.** The reduced-form VAR is
assumed invariant to the policy change being contemplated, which is precisely
what the Lucas critique denies. The exercise is informative about the model's
internal propagation, not about what the economy would truly have done, and
the further the counterfactual departs from historical experience the weaker
the claim. `C.ratio` — `||delta||` against a typical draw — is a usable gauge.

## Uncertainty bands

The paper reports conditional means only; Waggoner and Zha (1999) is largely
about the distribution around them. `wzBands` implements it.

**Shock uncertainty.** With `dy = Mstr eps`, `eps ~ N(0, I)` and the scenario
as `R eps = r`, conditioning a Gaussian on an affine restriction gives
`eps | R eps = r ~ N(pinv(R) r, I - pinv(R) R)`. The mean is the point
solution already computed; the covariance is the orthogonal projector onto
`null(R)`, applied as `z - pinv(R)(R z)` through the same SVD factors rather
than ever being formed. Because `R (I - pinv(R) R) = 0`, every draw satisfies
the constraint exactly and the constrained variable's band has zero width.

**Parameter uncertainty.** A residual bootstrap, with `R` and `r` rebuilt
inside the loop: each draw moves `A`, `Sigma_u`, `D`, the MA coefficients and
the baseline, so the constraint must be re-imposed rather than the point
solution reused. The DGP is bias-adjusted following Kilian (1998) — the bias
is estimated by a first pass and subtracted, shrunk if that makes the
companion explosive. When pandemic dummies are active, residuals are resampled
from the non-dummied months only, so the 2020 observations do not come back in
through the bootstrap.

### The policy rate series

The default is `IR3TIB01CAM156N`, the OECD 3-month interbank rate, not the
Bank Rate `IRSTCB01CAM156N` the paper uses. The Bank Rate stopped updating:
FRED's migration to the new OECD Data Explorer carried only part of the Main
Economic Indicators across, and the test is the Notes block on the series page
— `OECD Data Filters: REF_AREA: ...` means migrated and live, the older
`OECD Descriptor ID:` style means frozen.

This is a concept change, not just a longer series. The Bank Rate is a policy
instrument moving in discrete 25bp steps; the interbank rate is
market-determined and embeds expectations of policy over the coming quarter.
It should help the row where the replication was furthest off — monetary
policy explaining 86–90% of its own forecast-error variance at long horizons
against the paper's 71% and 62% — but it changes what the tight-policy
scenario means, from a policy action to a statement about market pricing.
Revert with `opts.ids.rate = 'IRSTCB01CAM156N'`.

**Reading the fan chart — one trap.** The point baseline computed from the OLS
estimate is *not* the centre of the bootstrap distribution. The DGP is
bias-adjusted and explosive replications are rejected, so the baseline itself
shifts across draws; on this data the bootstrap median baseline sits visibly
away from the OLS baseline for the persistent variables (unemployment, the
policy rate) and essentially on top of it for the transitory ones. A gap
between the black baseline and the red scenario median therefore mixes the
scenario effect with that shift, and can appear even for a scenario that does
nothing. `figureBands` draws both baselines so the shift is visible; judge the
scenario against the dashed one. `figureDiffBands` differences within each
draw, which removes the issue entirely, and is the figure to read.

**Reading the output.** `B.scen(k).cond` is the band around the conditional
forecast. `B.scen(k).diff` is the band around scenario-minus-baseline, and it
carries parameter uncertainty *only*: with the parameters fixed that
difference is deterministic. `diff` is the object to read when asking whether
a scenario is distinguishable from the baseline, because the shock uncertainty
common to both cancels.

**A caveat with this data.** Persistence is close to a unit root, so the bias
adjustment pushes the DGP toward the unit circle and a large share of
bootstrap replications comes out explosive and is redrawn. `wzBands` warns
above a 20% rejection rate. The bands are then conditional on stationarity and
understate parameter uncertainty at long horizons; `biasCorrect = false`
lowers the rejection rate at the cost of a downward-biased persistence
estimate.

Plotting goes through `figureBands`, which transforms every draw with
`dispSeries` and takes quantiles afterwards — the display transforms are
nonlinear, so quantiles-then-transform would give the wrong band.

## A trap in the paper itself

The published **Figure 2 legend contradicts the Section 5 text and the Figure 3
legend**. Section 5 and Figure 3 both say purple = oil, orange = U.S. recession,
blue = unemployment, yellow = monetary policy; Figure 2's legend rotates three
of those four labels. (Check the BC Rate panel: the line that jumps to 6.5 and
settles at 4 is obviously the tight-MP path, but Figure 2's legend calls it
"U.S. Recession".) This package follows the text. If you compare curve by curve
against Figure 2's legend labels, three of four scenarios will look wrong when
they are not.

## Known gap: Figure 4

`figure4_mp_robustness.png` is **not** the paper's Figure 4 exercise. The paper
re-estimates on a pre-1992 sample, selects the 48 months of 1988–1991 with the
greatest policy-shock variability, and feeds those structural shocks to BK.
That needs data this package does not load. What is drawn instead takes the WZ
tightening innovations, maps them to structural space (`eps = D^{-1} u`), keeps
only the monetary-policy shock, and propagates via BK — the "tightening driven
purely by policy shocks" reading. To reproduce the paper, build the Appendix-C
`K x H` shock matrix and pass it to `bkConditional` as `Escen`.

## Files

| File | Role | Note equation |
|---|---|---|
| `estimateVAR.m` | OLS VAR(p) + companion form | (1)–(2) |
| `varMA.m` | reduced-form `Phi_j`, structural `Theta_j`, stacked maps | (5), (6′) |
| `baselineForecast.m` | unconditional forecast | (3)–(4) |
| `wzConditional.m` | Waggoner–Zha conditional forecast | (5)–(5″) |
| `bkConditional.m` | Baumeister–Kilian conditional forecast | (6)–(6′) |
| `varianceDecomp.m` | forecast-error variance decomposition | (7′)–(7″) |
| `checkTable1.m` | compare the decomposition with Table 1, gate the run | Appendix A |
| `seasonalityCheck.m` | month-dummy test for residual seasonality | — |
| `minnesotaDummies.m` | Minnesota prior as dummy observations | — |
| `varPosteriorDraw.m` | draw from the normal-inverse-Wishart posterior | — |
| `irfCholesky.m` | structural impulse responses, with bands | (6') |
| `figureIRF.m` | plot the responses, one figure per shock | — |
| `covidDummies.m` | pandemic-window dummies (diagnostic) | — |
| `wzBands.m` | shock + parameter uncertainty for the conditional forecasts | WZ (1999) |
| `wzCounterfactual.m` | historical counterfactual over an in-sample window | WZ (1999) |
| `figureBands.m` | fan chart per scenario | — |
| `figureDiffBands.m` | scenario minus baseline, differenced within draw | — |
| `selectLag.m` | AIC lag-order selection (diagnostic only) | Section 3 |
| `loadCanadaData.m` | assemble the 7 series, report spans | Section 3 |
| `fetchFRED.m`, `readFredCSV.m` | FRED CSV download / parse | — |
| `fetchStatCan.m` | download a StatCan table, filter it to one series | — |
| `readSeriesCSV.m` | reader for user-supplied series (long or StatCan wide) | — |
| `transformData.m` | log-differences / levels + level reconstruction | Section 3 |
| `buildScenarios.m` | constraint paths, spec-driven | Section 4 |
| `dispSeries.m`, `displayHelpers.m` | figure display conversions | Figures |
| `main_risk_scenarios.m` | driver + Figures 1–4 | — |
| `test_engines.m` | numerical validation harness | — |

## Method in one paragraph

`estimateVAR` fits the reduced form and builds the companion matrix `A`.
`baselineForecast` iterates `E[Y_{t+h}] = mu + A E[Y_{t+h-1}]`. `wzConditional`
stacks the deviation map `dy = M eps` (block-lower-triangular in
`Theta_j = Phi_j D`), imposes the scenario as `R eps = r`, and takes the
minimum-norm solution `pinv(R) r` via a rank-revealing SVD, then propagates it
to all variables. The norm is taken over the STRUCTURAL shocks — see below.
`bkConditional` instead takes a structural-shock path and propagates
`sum_j Theta_j eps`, with `D = chol(Sigma_u)`. `varianceDecomp` accumulates
`Theta_j Theta_j'` into per-shock FEV shares.

## Changelog

- Housing must now be a real monthly index; the spline proxy is opt-in and warns.
- `readSeriesCSV` reads StatCan's wide "Download as displayed" export directly —
  no pre-processing step, no `readtable`/`writetable` round trip.
- Terminal month of the panel is asserted; per-series spans are reported.
- Lag order fixed at `p = 3`; AIC kept as a diagnostic.
- Variance decomposition is computed and checked against Table 1 *before* any
  scenario, and gates the figures.
- U.S. IP and exchange-rate panels switched from annualised monthly growth to
  12-month growth.
- Monetary-policy scenario now pauses for exactly 8 months at the peak (was 9).
- `wzConditional` minimises the norm of the STRUCTURAL shocks rather than of
  the reduced-form innovations, which is what Doan-Litterman-Sims specifies and
  what "most likely shock combination" means. This changes every unconstrained
  variable in all four scenarios and in Figure 4. It uses a rank-revealing SVD
  and reports rank, conditioning and `||eps||`.
- `varianceDecomp` normalises via `bsxfun` (pre-R2016b and Octave safe).
- `readSeriesCSV` accepts a list of candidate row labels; `cpiFile` now routes
  a row label through (`cpiSeries`), which it previously did not.
- New `seasonalityCheck`, run before estimation.
- New `wzBands` / `figureBands`: shock and parameter uncertainty for the
  conditional forecasts, with the zero-width property of the constrained
  variable checked in `test_engines`.
- `rateFile` / `rateSeries` allow swapping the policy rate, mirroring the CPI
  and housing swaps. The FRED default `IRSTCB01CAM156N` is the Bank Rate, a
  step function that moves in discrete 25bp jumps and behaves close to
  exogenously; a market rate responds far more to macro developments.
- Fixed the figure x-axis: `addMonths` handed `datenum` an out-of-range month,
  which did not roll back a year, so the plots showed twelve months of history
  instead of forty-eight. The window is now `opts.figStart` (default 2019-01),
  matching the paper's Figures 1-3.
- `estimateVAR` and `selectLag` accept an exogenous block; new `covidDummies`
  supplies pandemic-window dummies as a diagnostic (`opts.covid`).
