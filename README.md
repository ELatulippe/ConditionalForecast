# Conditional Forecast — Canadian risk scenarios (MATLAB / Octave)

Replication of the Waggoner–Zha conditional-forecast methodology in
Moran, Stevanovic & Surprenant (2025, BoC SWP 2025-28): a monthly VAR for the
Canadian economy, a baseline forecast, risk scenarios, historical
counterfactuals, impulse responses and a forecast-error variance
decomposition. Runs on base MATLAB (R2014b+) and GNU Octave — no toolboxes.

## Folder layout

```
ConditionalForecast/
├── run_*.m                 # one runnable script per empirical exercise (top level)
├── run_all.m               # master driver: run every exercise, then build the PDF
├── freeze_panel.m          # one online run that builds the offline caches
├── build_report.m          # assemble figures + FEVD table into figures/report.{html,pdf}
├── scenario_config.m       # shared options struct for the scenario suite
├── scenarioOrigin.m        # shared "pass 1" origin read
├── paper/                  # the methods-and-tooling paper (paper.tex -> paper.pdf) + appendices
├── report.html             # static figure report (the "website" snapshot)
├── src/                    # all estimation / plotting helper functions
├── tests/                  # regression_test.m  (+ golden/ snapshot, created on first run)
├── data/                   # raw INPUT CSVs you supply (StatCan, BoC commodity index)
├── cache/                  # generated panel_*.mat caches (delete to rebuild)
├── figures/                # output figures + report.pdf (created on run)
└── docs/                   # the paper, the original README, Running_code.txt
```

`data/` holds only the inputs you provide; everything the code generates lands
in `cache/` and `figures/`. Keep the eight `run_*.m` scripts and the two
`scenario_*.m` helpers at the top level.

## Running

From MATLAB or Octave, `cd` to this folder and run whichever exercise you want:

```matlab
run_main               % baseline forecast + fan charts
run_single_scenarios   % four single-variable risk scenarios  -> Rb
run_joint_scenarios    % six multi-variable scenarios          -> Rj
run_temporary_shocks   % partial-horizon / transitory shocks   -> Rp
run_soft_conditions    % distributional (.sd) conditioning      -> Rs
run_counterfactuals    % in-sample counterfactuals (tighter BoC; no oil spike)
run_impulse_responses  % Cholesky IRFs (larger panel)          -> R12, I12
run_variance_decomp    % FEVD table + grouped area charts
```

Each script begins with

```matlab
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir); addpath(fullfile(thisDir,'src'));
```

so it puts `src/` on the path itself — you do **not** have to `addpath`
manually, and the scripts work from any current directory as long as they stay
at the project root. The scripts are independent: run them in any order.

To run everything and produce a single report of all the results:

```matlab
run_all       % runs every exercise, then calls build_report
```

`build_report` writes a self-contained `figures/report.html` (images embedded;
opens in any browser and prints to PDF with Ctrl+P — **no LaTeX needed**). If
your MATLAB is R2021b+ it also writes `figures/report.pdf` directly via
`exportgraphics`. A `report.tex` is always written too and is compiled to PDF
only if `pdflatex` happens to be on your PATH. So on a machine without LaTeX,
open `report.html`.

## Interactive dashboard

```matlab
build_dashboard        % writes dashboard/dashboard.html (self-contained)
```

`build_dashboard` estimates the VAR offline, computes the impulse responses and
their bands, embeds them (plus the scenario fan-chart PNGs already in
`figures/`) into a single self-contained `dashboard/dashboard.html`. Open it in
any browser — no server, no internet, no LaTeX. It has four tabs:

- **Impulse-response explorer** — pick the structural shock, drag its size
  (in standard-deviation multiples), set the horizon, and toggle
  level-cumulation and the uncertainty band. The responses are recomputed and
  redrawn *live in the browser* from the exported MA coefficients and Cholesky
  factor, so it is a true "what does a bigger oil shock do?" tool, not a static
  picture.
- **Scenario builder (live)** — impose your own Waggoner–Zha conditional
  forecast in the browser. Add one or more constraints (ramp / hold / pulse /
  policy / replay, with optional partial-horizon window and soft-condition sd),
  and the baseline-vs-scenario *level* paths for every variable recompute as you
  change the inputs. A feasibility banner reports the shock norm ‖ε‖ — how far
  the scenario sits in the model's tail. The whole conditional-forecast solve
  (build the stacked structural map from the MA coefficients, find the
  minimum-norm shocks satisfying the constraints, propagate, reconstruct levels)
  runs client-side and was cross-checked against MATLAB's `wzConditional` to
  about 1e-8. Tick **Show conditional band** to overlay the uncertainty around
  the scenario: it Monte-Carlos the conditional predictive distribution
  (`eps | constraint ~ N(êps, I − pinv(R)R)`) live, so the constrained variable
  has zero band width and the free variables fan out — the *shock* component of
  `wzBands`, computed without any bootstrap. (Estimation/parameter uncertainty
  is not included live — see the note below.)
- **Counterfactuals (live)** — the in-sample exercise: pick an origin month and
  window, override one variable's realized path (shift by Δ, hold at a value, or
  freeze a log-difference variable's level), and optionally confine the change
  to that variable's own structural shock. The browser computes the
  Waggoner–Zha historical counterfactual and plots actual vs counterfactual
  (cumulative % for growth variables, levels for the rest) with the end-of-window
  gap. Cross-checked against MATLAB `wzCounterfactual` to ~1e-8. A Lucas-critique
  caveat is shown because the estimated VAR is held fixed.
- **Scenario gallery** — flip through the precomputed scenario fan charts
  (populate them first by running the scenario scripts or `run_all`).

### A note on uncertainty bands

The scenario band computed live captures **shock (innovation) uncertainty with
the parameters held fixed** — the honest "given the constraint, how uncertain is
the rest of the economy?" band, and it needs no bootstrap because it is a
closed-form Gaussian. It does **not** include **estimation (parameter)
uncertainty**, which requires drawing VAR parameters and re-solving the
conditional forecast for each draw (that is what `wzBands.m` does, and what the
precomputed fan-chart PNGs in the gallery show — full bands including both
sources). Making the *full* band live for arbitrary user scenarios would mean
exporting a set of posterior/bootstrap parameter draws and re-running the solve
over them in the browser; feasible but data-heavy, so it is intentionally left
to the offline `wzBands` path.

`dashboard/dashboard_demo.html` is a ready-to-open copy built from **synthetic**
data, so you can see the interface immediately without running MATLAB. Requires
MATLAB R2016b+ (`jsonencode`).

## Offline mode (reproducible, no FRED calls)

Most series are pulled live from FRED, but each run also caches the assembled
panel to `cache/`. Once a matching cache exists, `main_risk_scenarios` reads it
and returns **without contacting FRED**, so runs are offline by default when the
caches are present (the two shipped caches, `panel_k9.mat` and `panel_k13.mat`,
already match the configs). To (re)build them from scratch in one online run:

```matlab
freeze_panel  % fetch everything once, write cache/panel_k9.mat and panel_k13.mat
```

For a *hard* guarantee — an error instead of a silent refetch if a cache is
missing or stale — set `ob.offline = true` in any config. In that mode the
loader never touches the network and tells you to run `freeze_panel` if it
cannot find a usable cache.

## Tests

```matlab
addpath('tests'); regression_test
```

`regression_test` runs the pipeline once, offline, from `cache/panel_k9.mat` and
checks algebraic invariants that must hold for any data — `Sigma_u` symmetric
PD, `D*D' = Sigma_u`, the impact IRF equals `D`, every FEVD row sums to 1, total
forecast-error variance is non-decreasing in the horizon, Waggoner–Zha hard
conditioning hits its constraint, and a counterfactual with `r_cf = r_actual`
reproduces history exactly. It also snapshots a handful of scalars to
`tests/golden/golden_regression.mat` on first run and compares against them
afterwards, so an unintended numerical change is caught. Delete that file to
re-baseline after an intended change. It errors on any failure, so it works as a
CI hook: `matlab -batch "addpath('tests'); regression_test"`.

## What changed in the data extraction (vs. the flat Running_code.zip)

Previously every CSV and `.mat` cache sat in one folder and was opened relative
to the current directory. Now that inputs live in `data/` and caches in
`cache/`, the loader resolves those locations itself:

* **`main_risk_scenarios`** finds the project root (via `src/projectRoot.m`) and,
  unless you override them, sets `opts.dataDir = <root>/data` and
  `opts.cacheDir = <root>/cache`. `figures/` is likewise written under the root.
* **`loadCanadaData`** gained two options, `dataDir` and `cacheDir`, and three
  small internal helpers (`resolveData`, `resolveCache`, `isAbsolutePath`):
  * a **bare** input filename such as `'1810020501-eng.csv'` is looked up in
    `dataDir` first (and used if found there), otherwise left as-is;
  * a **bare** cache name such as `'panel_k9.mat'` is placed under `cacheDir`
    for both load and save (the folder is created if needed);
  * an **absolute path**, or any path that already carries a folder, passes
    through untouched.

Because of that, the config files keep bare filenames — `scenario_config.m`
still says `'1810020501-eng.csv'` and `'panel_k9.mat'` — and they resolve to
`data/` and `cache/` automatically. The behaviour is fully backward compatible:
with no `data/`/`cache/` folders present, the old "everything relative to pwd"
path is used.

If your CSVs live somewhere else, point the loader at them explicitly, e.g.

```matlab
ob.dataDir  = '/path/to/my/csvs';
ob.cacheDir = '/path/to/my/caches';
```

## Inputs you must supply (already in `data/`)

| File | Series | Used for |
|---|---|---|
| `1810020501-eng.csv` | StatCan New Housing Price Index | `house` |
| `1810000601-eng.csv` | StatCan CPI (seasonally adjusted) | `cpi` |
| `1410028701-eng.csv` | StatCan unemployment rate | `unemp` |
| `BCPI_MONTHLY-sd-1972-01-01.csv` | BoC non-energy commodity index | `bcne` |

The remaining series are pulled live from FRED (`source = 'fred'`), so a run
needs network access. To run entirely offline from pre-downloaded FRED CSVs,
use `source = 'local'` with `opts.csvDir` (or drop them in `data/`).

## Notes

* `src/test_engines.m` validates the numerical engines independently:
  `addpath('src'); test_engines`.
* Delete a `cache/panel_*.mat` file to force that panel to be reassembled after
  you change its variable set.
* `docs/Running_code.txt` is the original command log these scripts were built
  from; `docs/README_original.md` is the package's original documentation.
