"""Fast regression tests for mechanistic and numerical invariants."""
import unittest

import numpy as np
from hld_model import (
    Params, from_aggregates, integrate, invariant_errors,
    progress_qssa, v_qssa,
)


class ModelTests(unittest.TestCase):
    def test_aggregates_recovered(self):
        for gamma in (1.0, 1.5, 3.0):
            with self.subTest(gamma=gamma):
                p = from_aggregates(6.78, 16.9, 7.54, gamma=gamma)
                # gamma=1 is only an approximate fast-downstream limit.
                rtol = 5e-5 if gamma == 1.0 else 1e-10
                self.assertTrue(np.isclose(p.Km, 6.78, rtol=rtol))
                self.assertTrue(np.isclose(p.kcat, 16.9, rtol=rtol))
                self.assertTrue(np.isclose(p.KSI, 7.54, rtol=rtol))
                if gamma > 1.0:
                    self.assertAlmostEqual(p.gamma, gamma, places=8)

    def test_mass_conservation_and_nearly_nonnegative_solution(self):
        for gamma in (1.0, 3.0):
            with self.subTest(gamma=gamma):
                p = from_aggregates(6.78, 16.9, 7.54, gamma=gamma)
                sol = integrate(p, 10.0, np.linspace(0, 2, 101))
                e_err, s_err = invariant_errors(p, sol, 10.0)
                self.assertLess(e_err, 1e-10)
                self.assertLess(s_err, 1e-9)
                self.assertGreaterEqual(float(np.min(sol.y)), -1e-9)

    def test_rate_law_zero_and_no_inhibition(self):
        s = np.array([0.0, 1.0, 5.0])
        np.testing.assert_allclose(v_qssa(s, 2.0, 3.0),
                                   2.0 * s / (3.0 + s))
        self.assertEqual(v_qssa(0.0, 2.0, 3.0, 7.0), 0.0)

    def test_progress_with_no_enzyme_or_spontaneous_hydrolysis(self):
        p = from_aggregates(6.78, 16.9, 7.54, E0=0.0)
        t, y = progress_qssa(p, 5.0, np.linspace(0.0, 1.0, 11))
        self.assertEqual(len(t), 11)
        np.testing.assert_allclose(y, 0.0, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
