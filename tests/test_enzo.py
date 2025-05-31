import os
import sys

import yt
from stubs.libyt_stub import create_libyt_stub

import yt_libyt


def test_enzo_collapse_test_noncosmological():
    # Create libyt stub and make it importable
    simulation = "enzo"
    test_data_path = os.path.join(os.path.dirname(__file__), "data", simulation, "DD0000/DD0000")
    code_param_list = {
        "code_params": ["HydroMethod", "MultiSpecies", "DualEnergyFormalism"],
        "method": (lambda ds, code_param: ds.parameters[code_param]),
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
    slc.save()

    # Compare to post-processing results
    ds_post = yt.load(test_data_path)
    slc_post = yt.SlicePlot(ds_post, "z", ("enzo", "Density"))
    slc_post.save()

    assert sum(slc.data_source["enzo", "Density"] - slc_post.data_source["enzo", "Density"]) == 0.0
