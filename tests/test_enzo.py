import os
import sys

from stubs.libyt_stub import create_libyt_stub

import yt_libyt


def test_enzo():
    # Create libyt stub and make it importable
    simulation = "enzo"
    test_data_path = os.path.join(os.path.dirname(__file__), "data", simulation, "DD0000/DD0000")
    code_param_list = []
    field_list = {}
    particle_list = {}

    libyt_stub = create_libyt_stub(
        simulation, test_data_path, code_param_list, field_list, particle_list
    )
    sys.modules["libyt"] = libyt_stub

    # Run yt_libyt code
    ds = yt_libyt.libytDataset()
    ds.print_stats()
