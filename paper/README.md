# Paper

A methods-and-tooling paper that moves beyond the replication: its subject is the
*methods* — Waggoner–Zha conditional forecasting with uncertainty, soft and joint
scenarios, and historical counterfactuals — as implemented in this project.

## Files

| File | What it is |
|---|---|
| `paper.tex` | The paper. Main body (concise, intuition + key equations) + two appendices. Compiles to `paper.pdf`. |
| `app_original.tex` | Appendix A — the original equation-by-equation derivations from the source paper, `\input` by `paper.tex`. Auto-generated from `derivation_companion.tex` (sections demoted, labels namespaced). |
| `derivation_companion.tex` | The standalone derivation companion you provided, kept verbatim for reference (it also compiles on its own). |
| `paper.pdf` | The compiled paper. |

## Structure of the paper

- **Main body** (methods & tooling): the forecasting environment, conditional
  forecasts as a shock restriction, building scenarios (single / partial-horizon
  / joint / soft), **uncertainty** (the differentiator — closed-form shock band +
  simulated parameter band), historical counterfactuals, identification and the
  unemployment/employment caveat, the reproducible pipeline and interactive
  dashboard, and out-of-sample evaluation as future work.
- **Appendix A** (`app_original.tex`): the source model's derivations kept intact.
- **Appendix B**: the new derivations — structural vs reduced-form norm, the
  scenario-to-units map, joint/partial-horizon constraints, soft conditioning
  (Gaussian update), the shock-uncertainty band (with the zero-width proof for the
  constrained variable), the parameter-uncertainty band (posterior/bootstrap,
  Kilian bias adjustment), and the historical counterfactual.

## Compiling

```bash
pdflatex paper.tex
pdflatex paper.tex     # second pass for the table of contents and cross-refs
```

No `bibtex` step is needed (the bibliography is inline). No packages beyond a
standard TeX distribution.

## Figures

Figures are pulled from the pipeline's output in `../figures/` (the same PNGs the
HTML report embeds). The `\projfig` macro uses the real PNG when it is present and
otherwise draws a labelled placeholder box, so the paper compiles even before you
have run the pipeline. To populate the figures, run the scenario scripts (or
`run_all`) so `figures/` is filled, then recompile. The interactive versions of
these figures live in the dashboard (`../dashboard/`), and the static report is
`../report.html`.

## Regenerating Appendix A

If you edit `derivation_companion.tex` and want Appendix A to follow, re-run the
one-liner that produced `app_original.tex` (it strips the companion's preamble,
demotes its `\section`s to `\subsection`s, and namespaces its equation labels so
they don't clash with the main text).
