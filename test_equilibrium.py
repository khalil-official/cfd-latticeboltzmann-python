"""Check physical moments of the optimized D2Q9 equilibrium kernel."""
import unittest

import matplotlib
matplotlib.use("Agg")
import numpy as np

from latticeboltzmann import equilibrium


class EquilibriumTests(unittest.TestCase):
    def setUp(self):
        self.cxs = np.array([0, 0, 1, 1, 1, 0, -1, -1, -1])
        self.cys = np.array([0, 1, 1, 0, -1, -1, -1, 0, 1])
        self.weights = np.array([4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36])
        rng = np.random.default_rng(12)
        self.rho = rng.uniform(0.5, 2.0, (11, 17))
        self.ux = rng.uniform(-0.1, 0.1, self.rho.shape)
        self.uy = rng.uniform(-0.1, 0.1, self.rho.shape)

    def compute(self, out):
        return equilibrium(self.rho, self.ux, self.uy,
                           self.cxs, self.cys, self.weights, out)

    def test_equilibrium_preserves_density_and_momentum(self):
        populations = self.compute(np.empty((*self.rho.shape, 9)))
        np.testing.assert_allclose(populations.sum(axis=2), self.rho, rtol=1e-14)
        np.testing.assert_allclose((populations * self.cxs).sum(axis=2),
                                   self.rho * self.ux, rtol=1e-13, atol=1e-15)
        np.testing.assert_allclose((populations * self.cys).sum(axis=2),
                                   self.rho * self.uy, rtol=1e-13, atol=1e-15)

    def test_rest_state_has_isotropic_pressure_tensor(self):
        self.ux.fill(0)
        self.uy.fill(0)
        populations = self.compute(np.empty((*self.rho.shape, 9)))
        np.testing.assert_allclose((populations * self.cxs**2).sum(axis=2), self.rho / 3)
        np.testing.assert_allclose((populations * self.cys**2).sum(axis=2), self.rho / 3)
        np.testing.assert_allclose((populations * self.cxs * self.cys).sum(axis=2),
                                   0.0, atol=1e-15)

    def test_reused_buffer_is_fully_overwritten(self):
        out = np.full((*self.rho.shape, 9), np.nan)
        self.assertIs(self.compute(out), out)
        self.assertTrue(np.isfinite(out).all())
        self.rho *= 2
        self.assertIs(self.compute(out), out)
        np.testing.assert_allclose(out.sum(axis=2), self.rho, rtol=1e-14)


if __name__ == "__main__":
    unittest.main()
