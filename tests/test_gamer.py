import os
import sys

from stubs.libyt_stub import create_libyt_stub

import yt_libyt


def test_gamer_plummer():
    # Create libyt stub and make it importable
    simulation = "gamer"
    test_data_path = os.path.join(os.path.dirname(__file__), "data", simulation, "Plummer_000000")
    code_param_list = {
        "code_params": ["mhd", "gamma", "mu", "srhd"],
        "method": (lambda ds, code_param: getattr(ds, code_param)),
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

    libyt_stub = create_libyt_stub(
        simulation, test_data_path, code_param_list, field_list, particle_list
    )
    sys.modules["libyt"] = libyt_stub

    # Run yt_libyt code
    ds = yt_libyt.libytDataset()
    ds.print_stats()
