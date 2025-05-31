import types

import numpy as np
import yt

LIBYT_VERSION = (0, 2, 0)


def create_libyt_stub(
    simulation: str,
    fig_base_name: str,
    test_data: str,
    get_code_params: dict,
    field_list: dict,
    particle_list: dict,
    simulation_field_to_yt_field,
) -> types.ModuleType:
    """
    Returns a stub module that mimics libyt with a specific simulation.

    :note: yt store the data structure in 3d, which is if the simulation is 2d, an additional dimension is added.

          hierarchy is stored in 3d as well in a continuous array with 0 as the first grid, if dimensionality is 2,
          the third dimension is set to 1.
          grid_data stores the grid id starting from the index offset.

          When extracting the data using yt, the data is always in 3d, so if the simulation is 2d, we need to extract it.

    :param simulation: simulation name, e.g., "gamer", "enzo", etc.
    :param fig_base_name: figure base name
    :param test_data: the absolute path to the test data.
    :param get_code_params: the code parameters defined in the simulation frontend, and how to get it.
                            {
                             "code_params": [(param_name, default_value), ...],
                             "method": (function to get the parameter)
                             "expected_error": Exception type that is expected to be raised if the parameter is not found
                            }
    :param field_list: libyt-v0.2 defined field list
    :param particle_list: libyt-v0.2 defined particle list
    :param simulation_field_to_yt_field: mapping field_list and particle_list name to yt field name
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
        "fig_basename": fig_base_name,
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
        "field_list": None,
        "particle_list": None,
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
        try:
            stub.param_user[param[0]] = get_code_params["method"](ds, param[0])
        except get_code_params["expected_error"]:
            stub.param_user[param[0]] = param[1]

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

    # Fill in grid_data
    for gid in range(stub.param_yt["num_grids"]):
        if gid + stub.param_yt["index_offset"] not in stub.grid_data:
            stub.grid_data[gid + stub.param_yt["index_offset"]] = {}
        for field in field_list.keys():
            # TODO: assume cell-centered fields
            ghost_cells = field_list[field]["ghost_cell"]
            allocate_dim = stub.hierarchy["grid_dimensions"][gid][
                : stub.param_yt["dimensionality"]
            ].copy()
            if field_list[field]["contiguous_in_x"]:
                allocate_dim = np.flip(allocate_dim)
            for d in range(2 * stub.param_yt["dimensionality"]):
                allocate_dim[int(d / 2)] += ghost_cells[d]
            stub.grid_data[gid + stub.param_yt["index_offset"]][field] = (
                np.ones(allocate_dim, dtype=np.float64) * np.nan
            )

            if field_list[field]["contiguous_in_x"]:
                if stub.param_yt["dimensionality"] == 3:
                    stub.grid_data[gid + stub.param_yt["index_offset"]][field][
                        ghost_cells[0] : allocate_dim[0] - ghost_cells[1],
                        ghost_cells[2] : allocate_dim[1] - ghost_cells[3],
                        ghost_cells[4] : allocate_dim[2] - ghost_cells[5],
                    ] = (
                        ds.index.grids[gid][simulation_field_to_yt_field[field]]
                        .swapaxes(0, 2)
                        .in_base("code")
                    )
                elif stub.param_yt["dimensionality"] == 2:
                    stub.grid_data[gid + stub.param_yt["index_offset"]][field][
                        ghost_cells[0] : allocate_dim[0] - ghost_cells[1],
                        ghost_cells[2] : allocate_dim[1] - ghost_cells[3],
                    ] = np.squeeze(
                        ds.index.grids[gid][simulation_field_to_yt_field[field]].in_base("code")
                    ).swapaxes(
                        0, 1
                    )
                else:
                    stub.grid_data[gid + stub.param_yt["index_offset"]][field][
                        ghost_cells[0] : allocate_dim[0] - ghost_cells[1]
                    ] = ds.index.grids[gid][simulation_field_to_yt_field[field]].in_base("code")

    return stub
