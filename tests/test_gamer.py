import os
import sys

import yt
from stubs.libyt_v_0_2_stub import create_libyt_stub

import yt_libyt


def test_gamer_plummer():
    # Create libyt stub and make it importable
    simulation = "gamer"
    problem = "Plummer"
    fig_base_name = f"{simulation}-{problem}"
    test_data_path = os.path.join(
        os.path.dirname(__file__), "data", simulation, f"{problem}/plummer_000000"
    )
    code_param_list = {
        "code_params": [
            ("mhd", None),
            ("gamma", None),
            ("mu", None),
            ("srhd", None),
            ("opt_unit", 0),
        ],
        "method": (lambda ds, code_param: getattr(ds, code_param)),
        "expected_error": AttributeError,
    }
    field_list = {
        "Dens": {
            "attribute": ["", [], None],
            "field_type": "cell-centered",
            "contiguous_in_x": True,
            "ghost_cell": [0, 0, 0, 0, 0, 0],
        }
    }
    particle_list = {}
    simulation_field_to_yt_field = {"Dens": ("gas", "density")}

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
    slc = yt.SlicePlot(ds, "z", ("gamer", "Dens"))
    slc.save(f"{fig_base_name}-libyt.png")

    # Compare to post-processing results
    ds_post = yt.load(test_data_path)
    slc_post = yt.SlicePlot(ds_post, "z", ("gamer", "Dens"))
    slc_post.save(f"{fig_base_name}-post.png")
    fig_size = slc_post.data_source["gamer", "Dens"].size

    assert (
        sum(slc.data_source["gamer", "Dens"] - slc_post.data_source["gamer", "Dens"]) / fig_size
        < 1e-2
    )
