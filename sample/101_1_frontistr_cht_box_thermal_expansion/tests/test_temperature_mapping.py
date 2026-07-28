import tempfile
import unittest
from pathlib import Path
import sys

import numpy as np
import yaml

PYTHON_DIR = Path(__file__).resolve().parents[1] / "python"
sys.path.insert(0, str(PYTHON_DIR))

from openfoam_temperature import (  # noqa: E402
    _read_of_scalar_field,
    _read_of_vector_field,
    align_cell_centers_to_node_mesh,
    interpolate_to_nodes,
)


class OpenFoamFieldReaderTests(unittest.TestCase):
    def test_compact_scalar_field(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "T"
            path.write_text(
                "internalField nonuniform List<scalar> 3(293.15 294.0 295.5);",
                encoding="utf-8",
            )
            np.testing.assert_allclose(_read_of_scalar_field(path), [293.15, 294.0, 295.5])

    def test_compact_vector_field(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "C"
            path.write_text(
                "internalField nonuniform List<vector> 2((0.4 0.4 0.1) (0.6 0.6 0.3));",
                encoding="utf-8",
            )
            np.testing.assert_allclose(
                _read_of_vector_field(path),
                [[0.4, 0.4, 0.1], [0.6, 0.6, 0.3]],
            )


class CoordinateMappingTests(unittest.TestCase):
    def test_openfoam_centered_box_is_shifted_to_frontistr_origin(self):
        cell_centers = np.array([[0.425, 0.425, 0.025], [0.575, 0.575, 0.375]])
        node_coords = np.array([[0.0, 0.0, 0.0], [0.2, 0.2, 0.4]])
        aligned, translation = align_cell_centers_to_node_mesh(cell_centers, node_coords)

        np.testing.assert_allclose(translation, [-0.4, -0.4, 0.0])
        np.testing.assert_allclose(aligned[0], [0.025, 0.025, 0.025])
        np.testing.assert_allclose(aligned[1], [0.175, 0.175, 0.375])

    def test_uniform_temperature_remains_uniform_after_interpolation(self):
        nodes = np.array([[0.0, 0.0, 0.0], [0.2, 0.2, 0.4]])
        centers = np.array([[0.025, 0.025, 0.025], [0.175, 0.175, 0.375]])
        temperatures = np.array([293.15, 293.15])

        np.testing.assert_allclose(
            interpolate_to_nodes(nodes, centers, temperatures, k=2),
            [293.15, 293.15],
        )

    def test_inverse_distance_weighted_temperature(self):
        nodes = np.array([[0.0, 0.0, 0.0]])
        centers = np.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        temperatures = np.array([300.0, 400.0])

        # weights=(1, 1/2), so T=(1*300 + 0.5*400)/(1+0.5)
        np.testing.assert_allclose(
            interpolate_to_nodes(nodes, centers, temperatures, k=2),
            [333.3333333333333],
        )


class MaterialConfigurationTests(unittest.TestCase):
    def test_material_yaml_contains_required_values(self):
        material_path = (
            Path(__file__).resolve().parents[1]
            / "config"
            / "material_properties_steel.yaml"
        )
        material = yaml.safe_load(material_path.read_text(encoding="utf-8"))

        required = {
            "young_modulus_Pa",
            "poisson_ratio",
            "density_kg_m3",
            "thermal_expansion_coeff_per_K",
            "reference_temperature_K",
        }
        self.assertTrue(required.issubset(material))
        for key in required:
            self.assertIsInstance(material[key], (int, float), key)
        self.assertEqual(material["reference_temperature_K"], 293.15)


if __name__ == "__main__":
    unittest.main()
