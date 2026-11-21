#!/usr/bin/env python3
"""
Reproduce every number and the figure of the accompanying analysis.

Usage:   python3 run_all.py            (writes results.txt and the figure)
         python3 run_all.py --quiet    (figure only)

Runtime: about 3 minutes on a laptop.  Random seeds are fixed, so the
tables reproduce exactly; the last digit of a bootstrap statistic may
differ across BLAS builds.
"""

import io
import subprocess
import sys
from contextlib import redirect_stdout

SECTIONS = [
    ("a1_qssa.py",
     "Invariants, QSSA validity, and the covalent-intermediate pool"),
    ("a2_design.py",
     "Substrate depletion within the initial-rate window"),
    ("a3_ident.py",
     "Parameter identifiability, reparametrisation, and design comparison"),
    ("fig_degeneracy.py",
     "Figure: the (Km, KSI) degeneracy and its removal"),
]


def main():
    quiet = "--quiet" in sys.argv
    chunks = []
    for script, label in SECTIONS:
        header = f"\n{'=' * 78}\n### {label}\n### source: {script}\n{'=' * 78}\n"
        if not quiet:
            print(header)
        proc = subprocess.run([sys.executable, script],
                              capture_output=True, text=True)
        if proc.returncode != 0:
            sys.stderr.write(f"\n{script} failed:\n{proc.stderr}\n")
            sys.exit(1)
        if not quiet:
            print(proc.stdout)
        chunks.append(header + proc.stdout)

    if not quiet:
        with open("results.txt", "w") as fh:
            fh.write("".join(chunks))
        print("\nWrote results.txt and fig_degeneracy.{pdf,png}")


if __name__ == "__main__":
    main()
