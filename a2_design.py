"""Analysis 2-3 -- depletion within the initial-rate window, and
   identifiability of (Vmax, Km, KSI) under the published design."""

import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import solve_ivp

rng = np.random.default_rng(20261103)

KM, KCAT, KSI = 6.78, 16.9, 7.54          # DphxA / COU2-Cl, Mion et al. 2026
GRID = np.array([0.8, 1.3, 2.1, 3.4, 5.5, 8.9,
                 14.3, 23.1, 37.3, 60.1, 80.0])   # log-spaced 0.8-80 uM


def law(S, V, Km, KSI):
    return V * S / (Km + S + S ** 2 / KSI)


# ======================================================================
# 2. Substrate depletion inside the measurement window
# ======================================================================
print("=" * 74)
print("2.  SUBSTRATE DEPLETION WITHIN THE INITIAL-RATE WINDOW")
print("=" * 74)
print("""
Fraction converted by time dt depends on Vmax and dt only through the
product Vmax*dt, so the table below is free of any assumption about E0.
Literature threshold for unbiased initial rates: 1% depletion.
""")


def frac_converted(S0, Vmax_dt):
    """Exact fraction converted, by integrating the reduced system."""
    def f(tau, y):            # tau = t/dt in [0,1], so dP/dtau = dt*v
        S = max(S0 - y[0], 0.0)
        return [Vmax_dt * S / (KM + S + S ** 2 / KSI)]
    sol = solve_ivp(f, (0, 1), [0.0], rtol=1e-10, atol=1e-14)
    return sol.y[0, -1] / S0


print(f"{'Vmax*dt (uM)':>13} |" + "".join(f"{s:>9.1f}" for s in GRID[:7]))
print("-" * 76)
for Vdt in [0.001, 0.01, 0.05, 0.1, 0.5, 1.0]:
    row = "".join(f"{100*frac_converted(s, Vdt):>8.2f}%" for s in GRID[:7])
    print(f"{Vdt:>13.3f} |{row}")

# largest Vmax*dt keeping depletion under threshold at the lowest point
for thr in [0.01, 0.05]:
    lo, hi = 1e-6, 10.0
    for _ in range(80):
        mid = np.sqrt(lo * hi)
        if frac_converted(GRID[0], mid) < thr:
            lo = mid
        else:
            hi = mid
    print(f"\n  To keep depletion below {100*thr:.0f}% at S0 = 0.8 uM:  "
          f"Vmax*dt < {lo:.4g} uM")
    for dt in [10, 30, 60, 300]:
        E0_max = lo / dt / KCAT
        print(f"      window {dt:4d} s  ->  Vmax < {lo/dt:.3e} uM/s  "
              f"->  E0 < {1e6*E0_max:.3f} pM")

# ======================================================================
# 3. Identifiability: what noise reproduces the published error bars?
# ======================================================================
print()
print("=" * 74)
print("3.  IDENTIFIABILITY UNDER THE PUBLISHED CONCENTRATION GRID")
print("=" * 74)
print("""
Published (DphxA / COU2-Cl):  Km 6.78 +/- 3.01 (44%),
kcat 16.9 +/- 5.10 (30%),  KSI 7.54 +/- 3.26 (43%).
Synthetic replicates: 11 concentrations x 3 replicates, proportional noise.
""")

VTRUE = 1.0          # Vmax is a scale factor; CVs and correlations are
                     # independent of its value
TRUE = np.array([VTRUE, KM, KSI])


def one_fit(sigma):
    S = np.repeat(GRID, 3)
    v = law(S, *TRUE)
    y = v * (1.0 + sigma * rng.standard_normal(S.size))
    try:
        popt, _ = curve_fit(law, S, y, p0=TRUE, maxfev=40000,
                            bounds=([0, 0, 0], [np.inf, np.inf, np.inf]))
        return popt
    except Exception:
        return None


print(f"{'noise':>7} {'CV(Vmax)':>10} {'CV(Km)':>9} {'CV(KSI)':>9} "
      f"{'r(Km,KSI)':>11} {'r(V,Km)':>9} {'r(V,KSI)':>10}")
print("-" * 74)
for sigma in [0.02, 0.05, 0.10, 0.15, 0.20, 0.30]:
    fits = [one_fit(sigma) for _ in range(600)]
    F = np.array([f for f in fits if f is not None])
    cv = F.std(axis=0, ddof=1) / F.mean(axis=0)
    C = np.corrcoef(F.T)
    print(f"{100*sigma:6.0f}% {100*cv[0]:9.1f}% {100*cv[1]:8.1f}% "
          f"{100*cv[2]:8.1f}% {C[1,2]:11.3f} {C[0,1]:9.3f} {C[0,2]:10.3f}")

print("""
Natural-scale correlations and CVs become unstable when fits run away;
these moments should not be interpreted as invariant across noise levels.
See a3_ident.py for bounded, log-parameter fits and a design comparison.
""")

# identifiable combinations, same simulations
print("Reparametrisation check at 15% noise "
      "(CV of derived quantities, 600 fits):")
F = np.array([f for f in (one_fit(0.15) for _ in range(600)) if f is not None])
V, Km, Ksi = F[:, 0], F[:, 1], F[:, 2]
derived = {
    "Vmax/Km   (low-S arm)": V / Km,
    "Vmax*KSI  (high-S arm)": V * Ksi,
    "S* = sqrt(Km*KSI)": np.sqrt(Km * Ksi),
    "Km                 ": Km,
    "KSI                ": Ksi,
}
for name, vals in derived.items():
    print(f"    {name}  CV = {100*vals.std(ddof=1)/vals.mean():6.1f}%")
