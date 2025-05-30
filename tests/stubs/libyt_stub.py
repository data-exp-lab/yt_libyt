import types

import numpy as np
import yt

LIBYT_VERSION = (0, 2, 0)


def create_libyt_stub(
    simulation: str, test_data: str, get_code_params: dict, field_list: dict, particle_list: dict
) -> types.ModuleType:
    """
    Returns a stub module that mimics libyt with a specific simulation.
    :param simulation: simulation name, e.g., "gamer", "enzo", etc.
    :param test_data: the absolute path to the test data.
    :param get_code_params: the code parameters defined in the simulation frontend, and how to get it.
    :param field_list: libyt-v0.2 defined field list
    :param particle_list: libyt-v0.2 defined particle list
    """
    # Mock libyt module based on libyt version 0.x.0
    stub = types.ModuleType("libyt")
    stub.libyt_info = {
        "version": LIBYT_VERSION,
        "SERIAL_MODE": True,
        "INTERACTIVE_MODE": False,
        "JUPYTER_KERNEL": False,
        "SUPPORT_TIMER": False,
    }
    stub.param_yt = {
        "frontend": simulation,
        "fig_basename": f"{simulation}0000",
        "current_time": None,
        "current_redshift": None,
        "omega_lambda": None,
        "omega_matter": None,
        "hubble_constant": None,
        "length_unit": None,
        "mass_unit": None,
        "time_unit": None,
        "magnetic_unit": None,
        "cosmological_simulation": None,
        "dimensionality": None,
        "refine_by": None,
        "velocity_unit": None,
        "domain_left_edge": None,
        "domain_right_edge": None,
        "periodicity": None,
        "domain_dimensions": None,
        "num_grids": None,
        "index_offset": None,
        "field_list": None,  # TODO
        "particle_list": None,  # TODO
    }
    stub.param_user = {}
    stub.hierarchy = {
        "grid_left_edge": None,
        "grid_right_edge": None,
        "grid_dimensions": None,
        "grid_parent_id": None,
        "grid_levels": None,
        "proc_num": None,
        # "par_count_list": None, # This is optional, which is a bad design choice.
    }
    stub.grid_data = {}
    stub.particle_data = {}

    ds = yt.load(test_data)
    # Fill in param_user
    for param in get_code_params["code_params"]:
        stub.param_user[param] = get_code_params["method"](ds, param)

    # Fill in param_yt
    for param in ds.__dict__.keys():
        if param in stub.param_yt and stub.param_yt[param] is None:
            stub.param_yt[param] = getattr(ds, param)
    if stub.param_yt["velocity_unit"] is None:
        stub.param_yt["velocity_unit"] = stub.param_yt["length_unit"] / stub.param_yt["time_unit"]
    stub.param_yt["domain_left_edge"] = ds.domain_left_edge
    stub.param_yt["domain_right_edge"] = ds.domain_right_edge
    stub.param_yt["periodicity"] = ds.periodicity
    stub.param_yt["domain_dimensions"] = ds.domain_dimensions
    stub.param_yt["num_grids"] = ds.index.num_grids
    stub.param_yt["index_offset"] = ds._index_class.grid._id_offset
    stub.param_yt["field_list"] = field_list
    stub.param_yt["particle_list"] = particle_list

    # Fill in hierarchy
    stub.hierarchy["grid_left_edge"] = ds.index.grid_left_edge
    stub.hierarchy["grid_right_edge"] = ds.index.grid_right_edge
    stub.hierarchy["grid_dimensions"] = ds.index.grid_dimensions
    stub.hierarchy["grid_levels"] = ds.index.grid_levels
    stub.hierarchy["proc_num"] = np.zeros((stub.param_yt["num_grids"], 1), dtype=np.int32)
    stub.hierarchy["grid_parent_id"] = np.ones((stub.param_yt["num_grids"], 1), dtype=np.int32) * -1
    for g in ds.index.grids:
        if g.Parent is not None:
            stub.hierarchy["grid_parent_id"][g.id - stub.param_yt["index_offset"]] = g.Parent.id

    return stub
