# Kinetic model and parameter identifiability for fluorogenic dehalogenase assays

Code accompanying the Supporting Information section *Kinetic model, assay design
and parameter identifiability* of:

> Mion, S., Mallet, M., Louis, Y., Fontanay, S., Ruaudel, F., Trujillo, L.,
> Bastian, F., Jacolot, M., Popowycz, F., Oger, P. M.
> *Coumarin-based blue excitable probe for functional screening of dehalogenase
> activity.*

This repository contains a mass-action model of a fluorogenic haloalkane
dehalogenase assay — including the covalent alkyl–enzyme intermediate, a
dead-end inhibitory complex and spontaneous hydrolysis of the probe — together
with the analyses built on it:

- how accurately the quasi-steady-state rate law describes the full system in
  the regime actually measured;
- how the fitted constants `Km`, `kcat` and `KSI` relate to the elementary rate
  constants of the catalytic mechanism;
- how much substrate is consumed inside the initial-rate measurement window;
- how well the concentration design determines each fitted parameter, and which
  combinations are better conditioned than the parameters themselves.

Everything runs on the aggregate kinetic parameters reported in the paper. No
experimental data files are required, and nothing here needs new measurements.

## Citing

Please cite both the article above and this software. `CITATION.cff` carries the
metadata; the archived release has its own DOI, which is the one to use for the
code.

## Requirements

Python 3.10 or later with NumPy, SciPy and Matplotlib:

```
pip install -r requirements.txt
```

The pinned versions in `requirements.txt` are the ones used to produce
`results.txt` and the figure. Any reasonably recent release should reproduce
them.

## Running

```
python3 run_all.py
```

Runs every analysis in order, writes `results.txt` with the full console output,
and produces `fig_degeneracy.pdf` and `fig_degeneracy.png`. Takes about three
minutes. Individual scripts can also be run on their own.

A reference copy of `results.txt` is committed, so a fresh run can be diffed
against it. Random seeds are fixed and the tables reproduce exactly; the last
digit of a bootstrap statistic may differ across BLAS builds.

## Which script produces what

| Script | Produces |
| --- | --- |
| `hld_model.py` | Module: the mass-action system, its invariants, the steady-state reduction, and the mapping between elementary and fitted constants. Imported by the others, not run directly. |
| `a1_qssa.py` | Invariant conservation, the comparison of full-system initial rates against the reduced law, the timescale separation, and the size of the covalent-intermediate pool. |
| `a2_design.py` | Converted fraction against concentration and against `Vmax*dt`, and the enzyme loadings that keep substrate depletion below threshold. Also contains a first-pass identifiability scan superseded by `a3_ident.py` (see below). |
| `a3_ident.py` | Identifiability of `(Vmax, Km, KSI)` under the concentration grid, the reparametrisation, and the comparison of alternative designs. |
| `fig_degeneracy.py` | The `(Km, KSI)` degeneracy figure, plus the geometric CVs of the identifiable combinations. |

## Implementation notes

**Integration.** `solve_ivp` with LSODA, `rtol=1e-10`, `atol=1e-14`. The two
conservation laws are checked along every trajectory and hold to a relative
error below `3e-15`; `hld_model.invariant_errors` computes them and is the
recommended solver diagnostic for schemes with sequestering complexes.

**Elementary constants.** These are not identifiable from steady-state data.
`from_aggregates` constructs a set consistent with given `Km`, `kcat` and
`KSI_app` for a chosen `gamma = 1 + k2/k3 + k2/k4`, with `k1` and `km5` fixed at
plausible values. The reported results do not depend on that choice: the
scenarios `gamma=1` and `gamma=3` give identical steady-state rates.

Note the distinction the module enforces. `Params.KSI` is the **apparent**
inhibition constant that appears in the fitted rate law, and `Params.KSI_true`
is the equilibrium dissociation constant of the dead-end complex. They differ by
`gamma`, because the dead-end complex competes against `ES` while `ES` holds
only a fraction `1/gamma` of the catalytically engaged enzyme.
`from_aggregates` takes the apparent one, since that is what a fit returns.

**Initial rates.** Extracted by least-squares slope over a window placed *after*
the binding transient (0.3–0.8 s). Including `t=0` mixes the transient into the
slope and inflates it; this matters when `gamma > 1`.

**Fitting.** `a3_ident.py` and `fig_degeneracy.py` fit on log-parameters with
proportional weighting, using `least_squares` with bounds. Log-space keeps the
parameters positive and makes the degenerate direction a straight line, which is
why the correlation estimates are stable. Fitting on the natural scale, as in
the first-pass scan inside `a2_design.py`, produces occasional runaway fits that
contaminate the moments — the two are kept side by side so the effect can be
compared.

**Statistics.** Dispersion is reported as a geometric coefficient of variation,
`exp(sd(log theta)) - 1`, comparable with a relative standard error and the
right scale for ratio parameters.

## License

MIT — see `LICENSE`.
