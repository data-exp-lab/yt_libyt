import os
import sys

import yt
from stubs.libyt_v_0_2_stub import create_libyt_stub

import yt_libyt


def test_enzo_isolated_galaxy():
    # Create libyt stub and make it importable
    simulation = "enzo"
    problem = "IsolatedGalaxy"
    fig_base_name = f"{simulation}-{problem}"
    test_data_path = os.path.join(
        os.path.dirname(__file__), "data", simulation, f"{problem}/galaxy0030/galaxy0030"
    )
    code_param_list = {
        "code_params": [
            ("HydroMethod", None),
            ("MultiSpecies", None),
            ("DualEnergyFormalism", None),
        ],
        "method": (lambda ds, code_param: ds.parameters[code_param]),
        "expected_error": KeyError,
    }
    field_list = {
        "Density": {
            "attribute": ["", [], None],
            "field_type": "cell-centered",
            "contiguous_in_x": True,
            "ghost_cell": [3, 3, 3, 3, 3, 3],
        }
    }
    particle_list = {}
    simulation_field_to_yt_field = {"Density": ("gas", "density")}

    libyt_stub = create_libyt_stub(
        simulation,
        fig_base_name,
        test_data_path,
        code_param_list,
        field_list,
        particle_list,
        simulation_field_to_yt_field,
    )
    sys.modules["libyt"] = libyt_stub

    # Run yt_libyt code
    ds = yt_libyt.libytDataset()
    slc = yt.SlicePlot(ds, "z", ("enzo", "Density"))
    slc.save(f"{fig_base_name}-libyt.png")

    # Compare to post-processing results
    ds_post = yt.load(test_data_path)
    slc_post = yt.SlicePlot(ds_post, "z", ("enzo", "Density"))
    slc_post.save(f"{fig_base_name}-post.png")
    fig_size = slc_post.data_source["enzo", "Density"].size

    assert (
        sum(slc.data_source["enzo", "Density"] - slc_post.data_source["enzo", "Density"]) / fig_size
        < 1e-2
    )


def test_enzo_kelvin_helmholtz():
    # Create libyt stub and make it importable
    simulation = "enzo"
    problem = "EnzoKelvinHelmholtz"
    fig_base_name = f"{simulation}-{problem}"
    test_data_path = os.path.join(
        os.path.dirname(__file__), "data", simulation, f"{problem}/DD0000/DD0000"
    )
    code_param_list = {
        "code_params": [
            ("HydroMethod", None),
            ("MultiSpecies", None),
            ("DualEnergyFormalism", None),
        ],
        "method": (lambda ds, code_param: ds.parameters[code_param]),
        "expected_error": KeyError,
    }
    field_list = {
        "Density": {
            "attribute": ["", [], None],
            "field_type": "cell-centered",
            "contiguous_in_x": True,
            "ghost_cell": [3, 3, 3, 3, 0, 0],
        }
    }
    particle_list = {}
    simulation_field_to_yt_field = {"Density": ("gas", "density")}

    libyt_stub = create_libyt_stub(
        simulation,
        fig_base_name,
        test_data_path,
        code_param_list,
        field_list,
        particle_list,
        simulation_field_to_yt_field,
    )
    sys.modules["libyt"] = libyt_stub

    # Run yt_libyt code
    ds = yt_libyt.libytDataset()
    slc = yt.SlicePlot(ds, "z", ("enzo", "Density"))
    slc.save(f"{fig_base_name}-libyt.png")

    # Compare to post-processing results
    ds_post = yt.load(test_data_path)
    slc_post = yt.SlicePlot(ds_post, "z", ("enzo", "Density"))
    slc_post.save(f"{fig_base_name}-post.png")
    fig_size = slc_post.data_source["enzo", "Density"].size

    assert (
        sum(slc.data_source["enzo", "Density"] - slc_post.data_source["enzo", "Density"]) / fig_size
        < 1e-2
    )
