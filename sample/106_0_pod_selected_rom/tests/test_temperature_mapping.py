"""IDW checks: analytic weights, exact matches, constants, linearity, axes."""
import unittest

import numpy as np

from dacore.temperature_mapping import build_idw_mapping, map_temperature


class TemperatureMappingTest(unittest.TestCase):
    def test_known_inverse_distance_weights(self):
        source = np.array([[1., 0, 0], [2., 0, 0], [4., 0, 0]])
        mapping = build_idw_mapping(source, [[0., 0, 0]], k=8)
        np.testing.assert_allclose(mapping.toarray(), [[4/7, 2/7, 1/7]])
        np.testing.assert_allclose(map_temperature(mapping, [25., 23., 21.]), [41.75/1.75])

    def test_exact_match_and_constant_preservation(self):
        source = [[0., 0, 0], [1., 0, 0], [2., 0, 0]]
        mapping = build_idw_mapping(source, [[1., 0, 0], [.2, 0, 0]])
        np.testing.assert_allclose(mapping.toarray()[0], [0, 1, 0])
        np.testing.assert_allclose(np.asarray(mapping.sum(axis=1)).ravel(), 1)
        np.testing.assert_allclose(map_temperature(mapping, np.full(3, 293.15)), 293.15)
        self.assertTrue((mapping.data >= 0).all())

    def test_linearity_and_time_axis(self):
        source = [[0., 0, 0], [1., 0, 0], [2., 0, 0]]
        mapping = build_idw_mapping(source, [[.2, 0, 0], [1.5, 0, 0]])
        modes = np.array([[1., 2.], [3., 4.], [5., 6.]])
        coefficients = np.array([.7, -.2])
        np.testing.assert_allclose(map_temperature(mapping, modes @ coefficients),
                                   map_temperature(mapping, modes) @ coefficients)
        np.testing.assert_allclose(map_temperature(mapping, modes.T, axis=1),
                                   map_temperature(mapping, modes).T)


if __name__ == "__main__":
    unittest.main()
