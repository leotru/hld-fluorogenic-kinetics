"""
Mass-action model of the COU2-Cl / COU4-Cl haloalkane dehalogenase assay.

Scheme (Mion et al. 2026, mechanism per Prokop et al. 2003):

    E + S  <--k1/km1-->  ES  --k2-->  EI  --k3-->  EP  --k4-->  E + P
    ES + S <--k5/km5-->  SES        (dead-end)
    S      --k0-->       P          (spontaneous hydrolysis)

State vector x = [S, ES, SES, EI, EP, P];  E from conservation.

Units: concentrations in uM, time in s.
"""

import numpy as np
from scipy.integrate import solve_ivp

# ----------------------------------------------------------------------
# Parameter container
# ----------------------------------------------------------------------


class Params:
    """Elementary rate constants, plus the aggregates they imply."""

    def __init__(self, k1, km1, k2, k3, k4, k5, km5, k0, E0):
        self.k1, self.km1 = k1, km1
        self.k2, self.k3, self.k4 = k2, k3, k4
        self.k5, self.km5 = k5, km5
        self.k0 = k0
        self.E0 = E0

    # --- aggregate (steady-state) parameters -------------------------
    @property
    def gamma(self):
        return 1.0 + self.k2 / self.k3 + self.k2 / self.k4

    @property
    def kcat(self):
        return self.k2 / self.gamma

    @property
    def Km(self):
        return (self.km1 + self.k2) / self.k1 * (self.kcat / self.k2)

    @property
    def Ks(self):
        """True dissociation constant of ES."""
        return self.km1 / self.k1

    @property
    def KSI_true(self):
        """Dissociation constant of the dead-end complex, km5/k5."""
        return np.inf if self.k5 == 0 else self.km5 / self.k5

    @property
    def KSI(self):
        """
        APPARENT inhibition constant appearing in the fitted rate law
            v = kcat*E0*S / (Km + S + S^2/KSI).

        Because the dead-end complex competes against ES only, while ES is
        a fraction 1/gamma of the catalytically engaged enzyme, the
        apparent constant is gamma times the true dissociation constant.
        """
        return self.gamma * self.KSI_true

    @property
    def Vmax(self):
        return self.kcat * self.E0

    def summary(self):
        return (f"kcat={self.kcat:.3f} /s  Km={self.Km:.3f} uM  "
                f"KSI_app={self.KSI:.3f} uM  KSI_true={self.KSI_true:.3f} uM  "
                f"Ks={self.Ks:.3f} uM  gamma={self.gamma:.3f}  "
                f"Km/Ks={self.Km/self.Ks:.3f}")


def from_aggregates(Km, kcat, KSI_app=np.inf, k1=10.0, km5=100.0, k0=0.0,
                    E0=1e-3, gamma=1.0):
    """
    Build elementary constants reproducing published aggregate parameters.

    Km, kcat and KSI_app are the quantities a Michaelis-Menten fit returns.
    gamma = 1 + k2/k3 + k2/k4 = k2/kcat controls how much covalent
    intermediate accumulates: gamma -> 1 means downstream steps are fast
    and none builds up; gamma = 3 with k2 = k3 = k4 means it does.

    The true dissociation constant of the dead-end complex is
    KSI_true = KSI_app / gamma, and the true dissociation constant of ES
    is Ks = km1/k1, which differs from Km likewise.
    """
    k2 = kcat * gamma
    if gamma == 1.0:
        k3 = k4 = 1e6                     # effectively instantaneous
    else:
        # k2/k3 + k2/k4 = gamma - 1, split evenly
        k3 = k4 = 2.0 * k2 / (gamma - 1.0)
    km1 = k1 * Km * (k2 / kcat) - k2
    if np.isinf(KSI_app):
        k5 = 0.0
    else:
        k5 = km5 / (KSI_app / gamma)      # km5/KSI_true
    return Params(k1, km1, k2, k3, k4, k5, km5, k0, E0)


# ----------------------------------------------------------------------
# Full mass-action system
# ----------------------------------------------------------------------


def rhs(t, x, p):
    S, ES, SES, EI, EP, P = x
    E = p.E0 - ES - SES - EI - EP
    bind = p.k1 * E * S
    unbind = p.km1 * ES
    inhib = p.k5 * S * ES
    deinhib = p.km5 * SES
    spont = p.k0 * S

    dS = -bind + unbind - inhib + deinhib - spont
    dES = bind - (p.km1 + p.k2) * ES - inhib + deinhib
    dSES = inhib - deinhib
    dEI = p.k2 * ES - p.k3 * EI
    dEP = p.k3 * EI - p.k4 * EP
    dP = p.k4 * EP + spont
    return [dS, dES, dSES, dEI, dEP, dP]


def integrate(p, S0, t_eval, rtol=1e-10, atol=1e-14):
    """Integrate the full system from t=0 with substrate added to free enzyme."""
    x0 = [S0, 0.0, 0.0, 0.0, 0.0, 0.0]
    sol = solve_ivp(rhs, (0.0, t_eval[-1]), x0, t_eval=t_eval,
                    args=(p,), method="LSODA", rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


def invariant_errors(p, sol, S0):
    """Relative violation of the two conservation laws along the trajectory."""
    S, ES, SES, EI, EP, P = sol.y
    E = p.E0 - ES - SES - EI - EP
    enz = E + ES + SES + EI + EP
    sub = S + ES + 2.0 * SES + EI + EP + P
    return (np.max(np.abs(enz - p.E0)) / p.E0,
            np.max(np.abs(sub - S0)) / S0)


# ----------------------------------------------------------------------
# Steady-state (QSSA) rate law and its progress curve
# ----------------------------------------------------------------------


def v_qssa(S, Vmax, Km, KSI=np.inf):
    """Michaelis-Menten with optional dead-end substrate inhibition."""
    denom = Km + S + (S ** 2 / KSI if np.isfinite(KSI) else 0.0)
    return Vmax * S / denom


def progress_qssa(p, S0, t_eval):
    """Integrate the reduced one-equation system dP/dt = v(S0-P) + k0*(S0-P)."""
    Vmax, Km, KSI, k0 = p.Vmax, p.Km, p.KSI, p.k0

    def f(t, y):
        S = max(S0 - y[0], 0.0)
        return [v_qssa(S, Vmax, Km, KSI) + k0 * S]

    sol = solve_ivp(f, (0.0, t_eval[-1]), [0.0], t_eval=t_eval,
                    method="LSODA", rtol=1e-10, atol=1e-14)
    return sol.t, sol.y[0]
