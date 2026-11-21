"""Figure: the (Km, KSI) degeneracy and its removal by reparametrisation."""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

rng = np.random.default_rng(20261103)

KM, KSI, VTRUE = 6.78, 7.54, 1.0
GRID = np.array([0.8, 1.3, 2.1, 3.4, 5.5, 8.9, 14.3, 23.1, 37.3, 60.1, 80.0])
LOGTRUE = np.log([VTRUE, KM, KSI])
SIGMA, NREP, NBOOT = 0.10, 3, 1500


def law(S, V, Km, Ksi):
    return V * S / (Km + S + S ** 2 / Ksi)


def fit_once():
    S = np.repeat(GRID, NREP)
    v = law(S, VTRUE, KM, KSI)
    y = v * (1.0 + SIGMA * rng.standard_normal(S.size))

    def resid(lp):
        V, Km, Ksi = np.exp(lp)
        return (law(S, V, Km, Ksi) - y) / v

    r = least_squares(resid, LOGTRUE, method="trf",
                      bounds=(np.log([1e-4, 1e-4, 1e-4]),
                              np.log([1e4, 1e5, 1e5])))
    return r.x if r.success else None


L = np.array([x for x in (fit_once() for _ in range(NBOOT)) if x is not None])
V, Km, Ksi = np.exp(L[:, 0]), np.exp(L[:, 1]), np.exp(L[:, 2])

# identifiable combinations, normalised to their true values
A = (V / Km) / (VTRUE / KM)                 # low-substrate arm
B = (V * Ksi) / (VTRUE * KSI)               # high-substrate arm

r_raw = np.corrcoef(np.log(Km), np.log(Ksi))[0, 1]
r_rep = np.corrcoef(np.log(A), np.log(B))[0, 1]

fig, ax = plt.subplots(1, 2, figsize=(9.2, 4.3))

# ---- panel A: the degenerate direction -----------------------------
ax[0].scatter(Km / KM, Ksi / KSI, s=6, alpha=0.30,
              color="#3b6ea5", edgecolors="none")
lim = (0.45, 2.3)
x = np.geomspace(*lim, 50)
ax[0].plot(x, 1.0 / x, color="#444444", lw=1.3, ls="--",
           label=r"constant $S^{*}=\sqrt{K_m K_{SI}}$")
ax[0].plot(1, 1, marker="+", ms=14, mew=2.2, color="#c0392b",
           linestyle="none", label="true value")
ax[0].set_xscale("log"); ax[0].set_yscale("log")
ax[0].set_xlim(*lim); ax[0].set_ylim(*lim)
ax[0].set_aspect("equal")
ax[0].set_xlabel(r"$K_m$  (relative to true value)")
ax[0].set_ylabel(r"$K_{SI}^{\mathrm{app}}$  (relative to true value)")
ax[0].set_title(f"(a)  fitted constants   $r = {r_raw:.3f}$", fontsize=10.5)
ax[0].legend(fontsize=8.0, loc="upper right", framealpha=0.92)

# ---- panel B: the identifiable combinations ------------------------
ax[1].scatter(A, B, s=6, alpha=0.30, color="#3b6ea5", edgecolors="none")
ax[1].plot(1, 1, marker="+", ms=14, mew=2.2, color="#c0392b",
           linestyle="none")
ax[1].set_xscale("log"); ax[1].set_yscale("log")
ax[1].set_xlim(*lim); ax[1].set_ylim(*lim)
ax[1].set_aspect("equal")
ax[1].set_xlabel(r"low-substrate arm $V_{\max}/K_m$  (relative)")
ax[1].set_ylabel(r"high-substrate arm $V_{\max}K_{SI}^{\mathrm{app}}$  (relative)")
ax[1].set_title(f"(b)  identifiable combinations   $r = {r_rep:.3f}$",
                fontsize=10.5)

for a in ax:
    a.grid(True, which="major", lw=0.4, alpha=0.35)
    a.tick_params(labelsize=9)

fig.tight_layout()
# metadata pinned so repeated runs produce byte-identical files
fig.savefig("fig_degeneracy.pdf", bbox_inches="tight",
            metadata={"CreationDate": None})
fig.savefig("fig_degeneracy.png", dpi=160, bbox_inches="tight",
            metadata={"Software": None})

print(f"n fits          : {len(L)}")
print(f"r(lnKm,lnKSI)   : {r_raw:+.4f}")
print(f"r(lnA,lnB)      : {r_rep:+.4f}")
for nm, q in [("Km", Km / KM), ("KSI", Ksi / KSI), ("Vmax/Km", A),
              ("Vmax*KSI", B), ("S*", np.sqrt(Km * Ksi) / np.sqrt(KM * KSI))]:
    print(f"  gCV({nm:9s}) = {100*np.expm1(np.std(np.log(q), ddof=1)):6.1f}%")
