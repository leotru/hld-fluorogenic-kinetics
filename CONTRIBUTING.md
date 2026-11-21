# Contributing

Contributions improving reproducibility, documentation, numerical stability,
and model validation are welcome.

1. Open an issue explaining the scientific or software change before large revisions.
2. Work from a feature branch, with a clear commit message.
3. Use Python 3.11+ and install `requirements.txt` in an isolated environment.
4. Run `python -m unittest discover -s tests -v` and, for results-changing
   changes, `python run_all.py`. Explain any differences to the tracked
   reference outputs in your pull request.
5. Describe assumptions, units, input provenance, random seeds and solver
   tolerances. Do not claim simulated data are experimental observations.
6. Do not add raw participant data, confidential research files, credentials,
   or materials you are not permitted to redistribute.

By contributing, you agree that your contribution can be distributed under
the repository's MIT license. Credit and scientific authorship should be
discussed explicitly when a contribution is substantive; a code contribution
alone does not automatically determine authorship of an associated article.
