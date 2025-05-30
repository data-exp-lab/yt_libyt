import os
import sys

from stubs.libyt_stub import create_libyt_stub


def test_gamer():
    # Create libyt stub
    test_data_path = os.path.join(os.path.dirname(__file__), "data", "gamer", "Plummer_000000")
    libyt_stub = create_libyt_stub("gamer", test_data_path)
    sys.modules["libyt"] = libyt_stub

    import libyt

    assert libyt.libyt_info["version"] == "0.2.0"
