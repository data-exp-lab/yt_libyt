import types

# import yt


def create_libyt_stub(simulation: str, test_data: str) -> types.ModuleType:
    """
    Returns a stub module that mimics libyt with a specific simulation.
    """
    stub = types.ModuleType("libyt")

    # ds = yt.load(test_data)

    # Mock dictionary
    stub.libyt_info = {
        "version": "0.2.0",
        "SERIAL_MODE": True,
        "INTERACTIVE_MODE": False,
        "JUPYTER_KERNEL": False,
        "SUPPORT_TIMER": False,
    }

    stub.param_yt = {}

    stub.param_user = {}

    stub.param_hierarchy = {}

    stub.grid_data = {}

    stub.particle_data = {}

    return stub
