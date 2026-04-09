import os
import numpy as np
import pytest
from nemaktis.lc_material import DirectorField, QTensorField


def field_equal_ref(field_vals, nref):
    if field_vals.shape[-1] == 3:
        return np.allclose(field_vals, nref)
    elif field_vals.shape[-1] == 6:
        nref = np.array(nref)
        Qref = ( 3*nref[None,:]*nref[:,None] - np.eye(3) ) / 2
        return np.allclose(field_vals, [Qref[0,0],Qref[1,1],Qref[2,2],Qref[0,1],Qref[0,2],Qref[1,2]])


def q_from_n(n):
    Q = 1.5*np.outer(n,n) - 0.5*np.eye(3)
    return np.array([Q[0,0], Q[1,1], Q[2,2], Q[0,1], Q[0,2], Q[1,2]])


class TestTensorField:
    mesh_dimensions = (9, 8, 7)
    mesh_lengths = (8.0, 7.0, 6.0)
    n0_ez = np.array([0.0, 0.0, 1.0])
    n0_tilted = np.array([1.0, 0.0, 1.0]) / np.sqrt(2)

    @pytest.fixture(params=[DirectorField, QTensorField])
    def empty_field(self, request):
        FieldType = request.param
        return FieldType(mesh_lengths=self.mesh_lengths, mesh_dimensions=self.mesh_dimensions)

    @pytest.fixture
    def ez_field(self, empty_field):
        empty_field.vals[...] = self.n0_ez if empty_field._Nv == 3 else q_from_n(self.n0_ez)
        return empty_field

    @pytest.fixture
    def tilted_field(self, empty_field):
        empty_field.vals[...] = self.n0_tilted if empty_field._Nv == 3 else q_from_n(self.n0_tilted)
        return empty_field

    # --- Constructor ---

    def test_constructor_mesh_params(self, empty_field):
        assert empty_field.get_mesh_lengths() == self.mesh_lengths
        assert empty_field.get_mesh_dimensions() == self.mesh_dimensions
        assert empty_field.vals.shape == self.mesh_dimensions[::-1] + (empty_field._Nv,)
        assert np.allclose(empty_field.vals, 0.0)

    @pytest.mark.parametrize("FieldType", [DirectorField, QTensorField])
    def test_constructor_wrong_params_raises(self, FieldType):
        with pytest.raises(Exception):
            FieldType(mesh_lengths=(1.0, 1.0, 1.0))

    # --- vals property ---

    def test_vals_setter_getter(self, ez_field):
        assert ez_field.vals.shape == self.mesh_dimensions[::-1] + (ez_field._Nv,)
        assert np.allclose(ez_field.vals[:, :, :, 2], 1.0)

    def test_vals_setter_wrong_shape_raises(self, empty_field):
        with pytest.raises(Exception):
            empty_field.vals = np.zeros((4, 4, 4, 7))  # wrong number of components

    # --- get_mesh_dimensions / lengths / spacings ---

    def test_get_mesh_dimensions(self, ez_field):
        assert ez_field.get_mesh_dimensions() == self.mesh_dimensions

    def test_get_mesh_lengths(self, ez_field):
        assert ez_field.get_mesh_lengths() == self.mesh_lengths

    def test_get_mesh_spacings(self, ez_field):
        assert np.allclose(
            np.array(ez_field.get_mesh_spacings()),
            np.array(self.mesh_lengths) / (np.array(self.mesh_dimensions) - 1))

    def test_get_pos(self, ez_field):
        assert np.allclose(ez_field.get_pos(0, 0, 0), - np.array(self.mesh_lengths) / 2)

    def test_get_n_vertices(self, ez_field):
        assert ez_field.get_n_vertices() == np.prod(self.mesh_dimensions)

    # --- init_from_func(s) ---

    def test_init_from_funcs(self, empty_field):
        funcs = [lambda x, y, z: x + y + z] * empty_field._Nv
        empty_field.init_from_funcs(*funcs)
        assert empty_field.vals.shape == self.mesh_dimensions[::-1] + (empty_field._Nv,)
        val = [func(*empty_field.get_pos(0, 0, 0)) for func in funcs]
        assert(np.allclose(empty_field.vals[0,0,0], val))

    def test_init_from_func(self, empty_field):
        def func(x, y, z):
            return np.repeat((x+y+z)[...,None], empty_field._Nv, axis=-1)
        empty_field.init_from_func(func)
        assert empty_field.vals.shape == self.mesh_dimensions[::-1] + (empty_field._Nv,)
        assert(np.allclose(empty_field.vals[0,0,0], func(*empty_field.get_pos(0, 0, 0))))

    # --- mask ---

    def test_set_mask_droplet(self, ez_field):
        ez_field.set_mask(mask_type="droplet")
        assert ez_field.mask_type == "droplet"
        assert ez_field.mask_vals.shape == self.mesh_dimensions[::-1]
        print(np.sum(ez_field.mask_vals>0),
            4/3*np.pi*(self.mesh_lengths[2]/2)**3 / np.prod(ez_field.get_mesh_spacings()))
        assert np.isclose(
            np.sum(ez_field.mask_vals>0) * np.prod(ez_field.get_mesh_spacings()),
            4/3*np.pi*(self.mesh_lengths[2]/2)**3, rtol=0.05)

    def test_set_mask_formula(self, ez_field):
        formula = str((self.mesh_lengths[2]/2)**2)+"-x**2-y**2-z**2"
        ez_field.set_mask(mask_type="formula", mask_formula=formula)
        assert ez_field.mask_type == "formula"
        assert ez_field.mask_formula == formula
        assert np.isclose(
            np.sum(ez_field.mask_vals>0) * np.prod(ez_field.get_mesh_spacings()),
            4/3*np.pi*(self.mesh_lengths[2]/2)**3, rtol=0.05)

    def test_set_mask_raw(self, ez_field):
        mask = np.random.rand(*self.mesh_dimensions[::-1])
        ez_field.set_mask(mask_type="raw", mask_ndarray=mask)
        assert ez_field.mask_type == "raw"
        assert np.allclose(ez_field.mask_vals, mask)

    def test_set_mask_invalid_type_raises(self, ez_field):
        with pytest.raises(Exception):
            ez_field.set_mask(mask_type="invalid_type")

    def test_delete_mask(self, ez_field):
        ez_field.set_mask(mask_type="droplet")
        ez_field.delete_mask()
        assert ez_field.mask_type is None
        assert ez_field.mask_vals is None
        assert ez_field.mask_formula is None

    # --- extend ---

    def test_extend_increases_mesh(self, ez_field):
        ez_field.extend(2.0, 2.0)
        new_dims = ez_field.get_mesh_dimensions()
        assert np.all(new_dims[:2] > self.mesh_dimensions[:2])
        assert new_dims[2] == self.mesh_dimensions[2]

    def test_extend_vals_shape(self, ez_field):
        ez_field.extend(2.0, 2.0)
        assert ez_field.vals.shape == ez_field.get_mesh_dimensions()[::-1] + (ez_field._Nv,)

    def test_extend_with_mask(self, ez_field):
        ez_field.set_mask(mask_type="droplet")
        ez_field.extend(2.0, 2.0)
        assert ez_field.mask_vals.shape == ez_field.get_mesh_dimensions()[::-1]

    # --- rotate_90deg ---

    @pytest.mark.parametrize("axis", ["x", "y", "z"])
    def test_rotate_90deg_vals_shape(self, ez_field, axis):
        ez_field.rotate_90deg(axis)
        assert ez_field.vals.shape == ez_field.get_mesh_dimensions()[::-1] + (ez_field._Nv,)

    def test_rotate_90deg_x(self, ez_field):
        print(ez_field.get_mesh_dimensions(), ez_field.vals.shape)    
        ez_field.rotate_90deg("x")
        print(ez_field.get_mesh_dimensions(), ez_field.vals.shape)       
        assert np.all(ez_field.get_mesh_dimensions() == np.array(self.mesh_dimensions)[[0, 2, 1]])    
        assert field_equal_ref(ez_field.vals, [0.0, -1.0, 0.0])

    def test_rotate_90deg_y(self, ez_field):
        ez_field.rotate_90deg("y")        
        assert np.all(ez_field.get_mesh_dimensions() == np.array(self.mesh_dimensions)[[2, 1, 0]])     
        assert field_equal_ref(ez_field.vals, [1.0, 0.0, 0.0])

    def test_rotate_90deg_z(self, ez_field):
        ez_field.rotate_90deg("z")
        assert np.all(ez_field.get_mesh_dimensions() == np.array(self.mesh_dimensions)[[1, 0, 2]])
        assert field_equal_ref(ez_field.vals, [0.0, 0.0, 1.0])

    def test_rotate_90deg_invalid_axis_raises(self, ez_field):
        with pytest.raises(Exception):
            ez_field.rotate_90deg("w")

    # --- rotate_180deg ---

    @pytest.mark.parametrize("axis", ["x", "y", "z"])
    def test_rotate_180deg_vals_shape(self, ez_field, axis):
        ez_field.rotate_180deg(axis)
        assert self.mesh_dimensions == ez_field.get_mesh_dimensions()
        assert ez_field.vals.shape == self.mesh_dimensions[::-1] + (ez_field._Nv,)

    @pytest.mark.parametrize("axis", ["x", "y"])
    def test_rotate_180deg_xy(self, ez_field, axis):
        ez_field.rotate_180deg(axis)
        assert field_equal_ref(ez_field.vals, [0.0, 0.0, -1.0])
    
    def test_rotate_180deg_z(self, ez_field):
        ez_field.rotate_180deg("z")
        assert field_equal_ref(ez_field.vals, [0.0, 0.0, 1.0])

    def test_rotate_180deg_invalid_axis_raises(self, ez_field):
        with pytest.raises(Exception):
            ez_field.rotate_180deg("w")

    # --- rotate ---

    @pytest.mark.parametrize("axis", ["x", "y", "z"])
    def test_rotate_preserves_shape(self, ez_field, axis):
        ez_field.rotate(axis, 45.0)
        assert ez_field.get_mesh_dimensions() == self.mesh_dimensions
        assert ez_field.vals.shape == self.mesh_dimensions[::-1] + (ez_field._Nv,)

    def test_rotate_valid(self, tilted_field):
        tilted_field.rotate("y", 45.0)
        assert field_equal_ref(tilted_field.vals, [1.0, 0.0, 0.0])

    def test_rotate_invalid_axis_raises(self, ez_field):
        with pytest.raises(Exception):
            ez_field.rotate("w", 45.0)

    # --- rescale_mesh ---

    def test_rescale_mesh_lengths(self, ez_field):
        ez_field.rescale_mesh(2)
        assert np.allclose(ez_field.get_mesh_lengths(), np.array(self.mesh_lengths) * 2)

    def test_rescale_mesh_spacings(self, ez_field):
        old_mesh_spacings = ez_field.get_mesh_spacings()
        ez_field.rescale_mesh(2)
        assert np.allclose(ez_field.get_mesh_spacings(), np.array(old_mesh_spacings) * 2)

    # --- save_to_vti ---

    def test_save_to_vti(self, ez_field, tmp_path):
        file_path = str(tmp_path / "test_field")
        ez_field.save_to_vti(file_path)
        assert os.path.exists(file_path + ".vti")

    def test_save_to_vti_with_extension(self, ez_field, tmp_path):
        file_path = str(tmp_path / "test_field.vti")
        ez_field.save_to_vti(file_path)
        assert os.path.exists(file_path)

    # --- vti_file constructor ---

    def test_constructor_from_vti(self, ez_field, tmp_path):
        file_path = str(tmp_path / "test_field")
        ez_field.save_to_vti(file_path)
        loaded = ez_field.__class__(vti_file=file_path + ".vti")
        assert loaded.get_mesh_dimensions() == ez_field.get_mesh_dimensions()
        assert np.allclose(loaded.vals, ez_field.vals)

    # --- constraints and conversion ---

    def test_tensor_constraints(self, empty_field):
        if empty_field._Nv == 3:
            empty_field.normalize() # should not raise for zero field
            empty_field.vals[...] = [3.0, 4.0, 5.0]
            empty_field.normalize()
            assert np.allclose(np.linalg.norm(empty_field.vals, axis=3), 1.0)
        else:
            empty_field.vals[...] = [1.0, 2.0, 3.0, 0.5, 0.5, 0.5]
            empty_field.apply_traceless_constraint()
            assert np.allclose(np.sum(empty_field.vals[..., :3], axis=-1), 0.0)

    def test_tensor_conversion(self, tilted_field):
        if tilted_field._Nv == 3:
            qfield = tilted_field.get_qtensor_field()
            assert np.allclose(qfield.vals, q_from_n(self.n0_tilted))
        else:
            nfield = tilted_field.get_director_field()
            nfield.vals *= np.sign(nfield.vals[...,[2]]*self.n0_tilted[2])
            assert np.allclose(nfield.vals, self.n0_tilted)
