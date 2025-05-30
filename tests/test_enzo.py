import os
import sys

from stubs.libyt_stub import create_libyt_stub


def test_enzo():
    # Create libyt stub
    simulation = "enzo"
    test_data_path = os.path.join(os.path.dirname(__file__), "data", simulation, "DD0000/DD0000")
    code_param_list = []

    libyt_stub = create_libyt_stub(simulation, test_data_path, code_param_list)
    sys.modules["libyt"] = libyt_stub

    import libyt

    assert libyt.libyt_info["version"] == "0.2.0"
    for key in libyt.param_yt:
        assert libyt.param_yt[key] is not None, f"Parameter {key} should not be None"
