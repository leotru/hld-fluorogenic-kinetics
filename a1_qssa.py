"""Analysis 1 -- invariants and QSSA validity for DphxA / COU2-Cl."""

import numpy as np
from hld_model import from_aggregates, integrate, invariant_errors, progress_qssa

# Published aggregates, Mion et al. 2026 Table 2, DphxA with COU2-Cl
KM, KCAT, KSI = 6.78, 16.9, 7.54

print("=" * 72)
print("DphxA / COU2-Cl   Km = 6.78 uM   kcat = 16.9 /s   KSI = 7.54 uM")
print("=" * 72)

for gamma, label in [(1.0, "gamma=1  (no intermediate accumulation)"),
                     (3.0, "gamma=3  (k2=k3=k4, intermediate accumulates)")]:
    p = from_aggregates(KM, KCAT, KSI, gamma=gamma, E0=1e-3)
    print(f"\n--- {label}")
    print("   ", p.summary())
    print(f"    elementary: k1={p.k1:.3g} km1={p.km1:.4g} k2={p.k2:.4g} "
          f"k3={p.k3:.4g} k4={p.k4:.4g}")

    # ---- invariants ------------------------------------------------
    t = np.linspace(0, 60, 601)
    sol = integrate(p, 10.0, t)
    e_err, s_err = invariant_errors(p, sol, 10.0)
    print(f"    invariants: enzyme {e_err:.2e}, substrate {s_err:.2e} "
          "(relative, max over trajectory)")

    # ---- QSSA vs full system, initial rates ------------------------
    print(f"\n    {'S0 (uM)':>9} {'v_full':>12} {'v_qssa':>12} {'rel.err':>9}")
    for S0 in [0.8, 2.0, 5.0, 10.0, 20.0, 40.0, 80.0]:
        # initial rate from the full system over a short window
        tt = np.linspace(0.3, 0.8, 51)
        s = integrate(p, S0, tt)
        P = s.y[5]
        # least-squares slope over the window
        v_full = np.polyfit(tt, P, 1)[0]
        v_red = p.Vmax * S0 / (p.Km + S0 + S0 ** 2 / p.KSI)
        print(f"    {S0:9.1f} {v_full:12.5e} {v_red:12.5e} "
              f"{100*(v_full-v_red)/v_red:8.2f}%")

    # ---- timescale separation --------------------------------------
    S0 = 10.0
    tau_fast = 1.0 / (p.k1 * S0 + p.km1 + p.k2 + p.k5 * S0 + p.km5)
    v = p.Vmax * S0 / (p.Km + S0 + S0 ** 2 / p.KSI)
    tau_depl = S0 / v
    print(f"\n    timescales at S0=10 uM: fast {tau_fast:.3e} s, "
          f"depletion {tau_depl:.3e} s, ratio {tau_depl/tau_fast:.3e}")

    # ---- how much intermediate sits in the steady state -------------
    tt = np.linspace(0, 5, 501)
    s = integrate(p, 10.0, tt)
    S_, ES_, SES_, EI_, EP_, P_ = s.y
    D = EI_ + EP_ + P_
    idx = len(tt) // 2
    print(f"    at t={tt[idx]:.2f}s:  P={P_[idx]:.5f}  D={D[idx]:.5f}  "
          f"D-P={D[idx]-P_[idx]:.3e} uM  "
          f"(= {100*(D[idx]-P_[idx])/p.E0:.1f}% of E0)")
    print(f"    predicted pool  E0*(1-kcat/k2) = "
          f"{p.E0*(1-p.kcat/p.k2):.3e} uM")
