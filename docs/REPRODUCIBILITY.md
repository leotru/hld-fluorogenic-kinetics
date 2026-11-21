# Reproducibility statement

## Research object and limitations

This is a computational software release that accompanies *Kinetic model,
assay design and parameter identifiability* in the Supporting Information
of the Mion et al. preprint (2026), DOI: 10.64898/2026.10.02.756220.
It contains Python source, reference console output, and a reference figure.
It **does not** contain experimental raw data or provide an independent
reconstruction of the preprint's complete empirical analysis.

The kinetic scenarios use aggregate published parameter values for
DphxA / COU2-Cl as listed in the source scripts. Elementary rate constants
are constructed to match aggregate constants but are **not inferred
uniquely** from experimental observations. The Monte Carlo studies use
fixed pseudorandom seeds, synthetic proportional perturbations and fitting
procedures specified in the source files.

## Procedure

From the repository root, with Python 3.11+:

```sh
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python run_all.py
```

`run_all.py` runs the four numerical analyses and overwrites
`results.txt`, `fig_degeneracy.pdf` and `fig_degeneracy.png`. If you need
to preserve the committed reference files, first work in a disposable
clone or check the diff before restoring them.

## Technical environment and interpretation

Pinned environment: NumPy 2.4.4, SciPy 1.17.1, Matplotlib 3.10.8;
supported Python major/minor versions should be verified via CI.
Integrator: `scipy.integrate.solve_ivp` with LSODA and explicitly defined
tolerances; fits use SciPy least squares. Numerical output can change
slightly between platforms because of floating-point evaluation and
numerical-library implementations. Images are regenerated, so compare
scientific quantities rather than expecting byte-identical PDF/PNG files.

## Verification criteria

- The CI suite checks kinetic aggregate consistency, conservation of total
  enzyme and substrate, approximate numerical positivity, the limiting
  rate law and a zero-enzyme reduced trajectory.
- The full workflow is exercised manually or on a published release.
- Inspect output errors and compare statistics with the committed reference
  results using justified tolerances; no unconditional byte-for-byte equality
  is claimed.

## Scientific distinction: practical vs structural identifiability

Strong fitted-parameter correlations in finite, noisy concentration designs
demonstrate **practical ill-conditioning**, not by themselves structural
non-identifiability of the rate-law parameters. The analyses here examine
local and Monte Carlo sensitivity for the chosen designs. Neither code nor
metadata should describe this as a general formal identifiability theorem.

## Permanent archiving

A public GitHub release can be archived automatically by Zenodo. After
Zenodo issues a real DOI, add the version-specific DOI to the metadata and
README in a subsequent commit. Do not insert guessed DOI values into a
released snapshot. A GitHub commit hash and a Zenodo version DOI together
identify the executable source used for a result.
