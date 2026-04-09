import numpy as np
import pytest

from nemaktis.light_propagator import OpticalFields


class TestOpticalFields:
    @pytest.fixture
    def optical_fields_params(self):
        return {
            "wavelengths": [0.55, 0.65],
            "max_NA_objective": 0.8,
            "max_NA_condenser": 0.4,
            "N_radial_wavevectors": 3,
            "mesh_lengths": (6.0, 4.0),
            "mesh_dimensions": (7, 5),
        }

    @pytest.fixture
    def make_optical_fields(self, optical_fields_params):
        def _make(**kwargs):
            params = dict(optical_fields_params)
            params.update(kwargs)
            return OpticalFields(**params)

        return _make

    @pytest.fixture
    def optical_fields(self, make_optical_fields):
        return make_optical_fields()

    def test_constructor_sets_basic_properties(self, optical_fields_params, optical_fields):
        assert np.allclose(optical_fields.get_wavelengths(), optical_fields_params["wavelengths"])
        assert optical_fields.get_mesh_dimensions() == optical_fields_params["mesh_dimensions"]
        assert optical_fields.get_mesh_lengths() == optical_fields_params["mesh_lengths"]
        assert optical_fields.z_focus == 0

    def test_constructor_builds_expected_wavevector_count(self, optical_fields_params, optical_fields):
        Nr = optical_fields_params["N_radial_wavevectors"]
        expected_nq = 1 + 3 * Nr * (Nr - 1)
        assert optical_fields.get_wavevectors().shape == (expected_nq, 2)
            
    def test_constructor_raises_with_wrong_mesh_dimensions_size(self, make_optical_fields):
        with pytest.raises(Exception):
            make_optical_fields(mesh_dimensions=(7, 5, 3))

    def test_constructor_raises_with_wrong_mesh_lengths_size(self, make_optical_fields):
        with pytest.raises(Exception):
            make_optical_fields(mesh_lengths=(6.0, 4.0, 2.0))

    def test_mesh_spacings_and_position_are_consistent(self, optical_fields):
        dx, dy = optical_fields.get_mesh_spacings()
        assert np.isclose(dx, 1.0)
        assert np.isclose(dy, 1.0)

        x0, y0 = optical_fields.get_pos(0, 0)
        assert np.isclose(x0, -3.0)
        assert np.isclose(y0, -2.0)

    def test_get_n_vertices_matches_dimensions(self, optical_fields):
        Nx, Ny = optical_fields.get_mesh_dimensions()
        assert optical_fields.get_n_vertices() == Nx * Ny

    def test_vals_setter_accepts_correct_shape(self, optical_fields):
        vals = np.ones(optical_fields.vals.shape, dtype=np.complex128)
        optical_fields.vals = vals
        assert np.allclose(optical_fields.vals, vals)

    def test_vals_setter_rejects_wrong_shape(self, optical_fields):
        wrong_shape = np.ones((1, 1, 4, 2, 2), dtype=np.complex128)
        with pytest.raises(Exception):
            optical_fields.vals = wrong_shape

    def test_copy_returns_independent_object(self, optical_fields):
        optical_fields.vals[...] = (1.0 + 2.0j)
        copy_fields = optical_fields.copy()

        assert copy_fields is not optical_fields
        assert np.allclose(copy_fields.vals, optical_fields.vals)

        copy_fields.vals[...] = 0.0
        assert not np.allclose(copy_fields.vals, optical_fields.vals)

    def test_update_NA_objective_clamps_to_bounds(self, optical_fields):
        optical_fields.update_NA_objective(-0.5)
        optical_fields.focus_fields(0)
        zero_na_vals = optical_fields.focused_vals.copy()

        optical_fields.update_NA_objective(10.0)
        optical_fields.focus_fields(0)
        max_na_vals = optical_fields.focused_vals

        assert zero_na_vals.shape == optical_fields.vals.shape
        assert max_na_vals.shape == optical_fields.vals.shape

    def test_focus_fields_updates_focus_and_produces_output(self, optical_fields):
        optical_fields.vals[...] = (1.0 + 0.0j)
        optical_fields.focus_fields(z_focus=2.5)

        assert np.isclose(optical_fields.z_focus, 2.5)
        assert optical_fields.focused_vals.shape == optical_fields.vals.shape

    def test_get_delta_qr_and_get_qr_index(self, optical_fields_params, optical_fields):
        expected_delta = (
            optical_fields_params["max_NA_condenser"]
            / (optical_fields_params["N_radial_wavevectors"] - 1)
        )
        assert np.isclose(optical_fields.get_delta_qr(), expected_delta)

        idx = optical_fields.get_qr_index(0.21)
        assert idx == int(np.ceil(0.21 / expected_delta))

    def test_koehler_1d_wavevectors_when_nx_is_one(self, make_optical_fields):
        fields = make_optical_fields(mesh_dimensions=(1, 5), koehler_1D=True, N_radial_wavevectors=4)
        assert fields.get_wavevectors().shape == (7, 2)

    def test_koehler_1d_wavevectors_when_ny_is_one(self, make_optical_fields):
        fields = make_optical_fields(mesh_dimensions=(5, 1), koehler_1D=True, N_radial_wavevectors=4)
        assert fields.get_wavevectors().shape == (7, 2)
