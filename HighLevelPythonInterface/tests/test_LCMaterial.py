import os
import numpy as np
import pytest
from nemaktis.lc_material import *


class TestLCMaterial:
    @pytest.fixture(params=[DirectorField, QTensorField])
    def make_lc_material(self, request):
        FieldType = request.param
        def _make(**kwargs):
            params = {
                "lc_field": FieldType(mesh_lengths=(6.0, 7.0, 8.0), mesh_dimensions=(7, 8, 9)),
                "ne": 1.7, "no": 1.5, "nhost": 1, "nin": 1, "nout": 1, "ne_imag": 0
            }
            params.update(kwargs)
            return LCMaterial(**params)

        return _make

    def test_constructor_required_parameters(self, make_lc_material):
        material = make_lc_material()

        assert isinstance(material.lc_field, DirectorField) or isinstance(material.lc_field, QTensorField)
        assert material.ne == 1.7
        assert material.no == 1.5

    def test_constructor_default_optional_parameters(self, make_lc_material):
        material = make_lc_material()

        assert material.nhost == 1
        assert material.nin == 1
        assert material.nout == 1
        assert material.ne_imag == 0
        assert material.iso_layer_indices == []
        assert material.iso_layer_thicknesses == []
        assert material._zfoc_NA_corr == 0

    def test_constructor_custom_optional_parameters(self, make_lc_material):
        material = make_lc_material(nhost=1.33, nin=1.45, nout=1.6, ne_imag=0.02)

        assert material.nhost == 1.33
        assert material.nin == 1.45
        assert material.nout == 1.6
        assert material.ne_imag == 0.02    

    def test_constructor_raises_with_wrong_lc_field(self, make_lc_material):
        with pytest.raises(Exception):
            make_lc_material(lc_field="not_a_tensor_field")

    def test_add_isotropic_layer_adds_one_layer(self, make_lc_material):
        material = make_lc_material()
        material.add_isotropic_layer(nlayer=1.52, thickness=4.0)
        assert material.iso_layer_indices == [1.52]
        assert material.iso_layer_thicknesses == [4.0]
        assert material._zfoc_NA_corr == 3*material.nout*4.0 / (14*1.52) * (1/material.nout**2 - 1/1.52**2)

    def test_add_isotropic_layer_adds_multiple_layers_in_order(self, make_lc_material):
        material = make_lc_material()
        material.add_isotropic_layer(nlayer=1.4, thickness=1.0)
        material.add_isotropic_layer(nlayer=1.5, thickness=2.0)
        material.add_isotropic_layer(nlayer=1.6, thickness=3.0)

        assert material.iso_layer_indices == [1.4, 1.5, 1.6]
        assert material.iso_layer_thicknesses == [1.0, 2.0, 3.0]
        assert material._zfoc_NA_corr != 0

    def test_layers_are_instance_specific(self, make_lc_material):
        material_1 = make_lc_material()
        material_2 = make_lc_material()

        material_1.add_isotropic_layer(nlayer=1.55, thickness=5.0)

        assert material_1.iso_layer_indices == [1.55]
        assert material_1.iso_layer_thicknesses == [5.0]
        assert material_2.iso_layer_indices == []
        assert material_2.iso_layer_thicknesses == []
