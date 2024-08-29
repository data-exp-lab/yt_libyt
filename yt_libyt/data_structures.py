"""
libyt-specific data structures



"""

# -----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

import importlib
import os
import time
import weakref
from collections import defaultdict

import numpy as np
import yt.frontends.api
from yt.data_objects.index_subobjects.grid_patch import AMRGridPatch
from yt.data_objects.static_output import Dataset
from yt.funcs import mylog, setdefaultattr
from yt.geometry.geometry_handler import YTDataChunk
from yt.geometry.grid_geometry_handler import GridIndex

from .fields import libytFieldInfo


class libytGrid(AMRGridPatch):
    # We will set _id_offset to libyt.param_yt["index_offset"]
    # when initializing libytHierarchy
    _id_offset = 0

    def __init__(self, id, index, level):
        AMRGridPatch.__init__(self, id, filename=None, index=index)

        # Might be redundant by setting these two, as they refer to the same thing.
        # MPI_rank : Refers to grid is currently belongs to which MPI rank.
        # filename : Refers to grid is currently belongs to which MPI rank. See line 101.
        self.MPI_rank = None
        self.Parent = None
        self.Children = []
        self.Level = level

    def __repr__(self):
        return "libytGrid_%09i (dimension = %s)" % (self.id, self.ActiveDimensions)


class libytHierarchy(GridIndex):
    grid = libytGrid
    libyt = None

    # _preload_implemented = True # Not sure about this option

    def __init__(self, ds, dataset_type="libyt"):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        self.directory = os.getcwd()
        self.float_type = "float64"
        self.libyt = self.dataset.libyt

        libytGrid._id_offset = self.libyt.param_yt["index_offset"]

        GridIndex.__init__(self, ds, dataset_type)

    def _detect_output_fields(self):
        if "field_list" in self.libyt.param_yt:
            field_list = self.libyt.param_yt["field_list"]
            self.field_list = [
                (self.libyt.param_yt["frontend"].lower(), v) for v in field_list.keys()
            ]
        else:
            mylog.debug("No field list \"libyt.param_yt['field_list']\".")

        if "particle_list" in self.libyt.param_yt:
            particle_list = self.libyt.param_yt["particle_list"]
            for ptype in particle_list.keys():
                attribute = particle_list[ptype]["attribute"]
                self.field_list += [(ptype, particle) for particle in attribute.keys()]
        else:
            mylog.debug("No particle list \"libyt.param_yt['particle_list']\".")

    def _count_grids(self):
        self.num_grids = self.libyt.param_yt["num_grids"]

    def _initialize_grid_arrays(self):
        # We don't need to initialize a new buffer for hierarchy,
        # since these array are initialized and set in libyt.
        pass

    def _parse_index(self):
        mylog.debug("Setting index array for %s grids.", self.num_grids)

        # get hierarchy information of all grids from libyt
        hierarchy = self.libyt.hierarchy
        self.grid_dimensions = hierarchy["grid_dimensions"]
        self.grid_left_edge = self.ds.arr(hierarchy["grid_left_edge"], "code_length")
        self.grid_right_edge = self.ds.arr(hierarchy["grid_right_edge"], "code_length")
        self.grid_levels = hierarchy["grid_levels"]

        # derive max_level
        self.max_level = self.grid_levels.max()

        # derive grid_particle_count from par_count_list.
        # par_count_list is created only if there is particles.
        if "par_count_list" in hierarchy:
            self.grid_particle_count = np.sum(hierarchy["par_count_list"], axis=1)
            self.grid_particle_count = self.grid_particle_count[..., np.newaxis]
        else:
            self.grid_particle_count = np.zeros((self.num_grids, 1), "int32")

        # Indicates which MPI rank it belongs to.
        self.proc_num = hierarchy["proc_num"]

        # allocate all grid objects
        self.grids = np.empty(self.num_grids, dtype="object")
        index_offset = self.libyt.param_yt["index_offset"]
        for index in range(self.num_grids):
            self.grids[index] = self.grid(index + index_offset, self, self.grid_levels.flat[index])
            self.grids[index].MPI_rank = self.proc_num[index, 0]
            self.grids[index].filename = "MPI_Rank_%07i" % (self.proc_num[index, 0])

    def _populate_grid_objects(self):
        # must flat it since 'grid_parent_id' has the dimension [num_grids][1]
        parent_list = self.libyt.hierarchy["grid_parent_id"].flat
        index_offset = self.libyt.param_yt["index_offset"]

        for index in range(self.num_grids):
            grid = self.grids[index]
            parent_id = parent_list[index]

            # set up the parent-children relationship if parent grid exist
            if parent_id >= index_offset:
                parent_grid = self.grids[parent_id - index_offset]
                grid.Parent = parent_grid
                parent_grid.Children.append(grid)

            # set up other grid attributes
            grid._prepare_grid()
            grid._setup_dx()

    def _chunk_io(self, dobj, cache=True, local_only=False):
        gfiles = defaultdict(list)
        gobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for g in gobjs:
            gfiles[g.filename].append(g)

        for fn in sorted(gfiles):
            if local_only:
                gobjs = [g for g in gfiles[fn] if g.MPI_rank == self.comm.rank]
                gfiles[fn] = gobjs
            gs = gfiles[fn]
            count = self._count_selection(dobj, gs)
            yield YTDataChunk(dobj, "io", gs, count, cache=cache)


class libytDataset(Dataset):
    _index_class = libytHierarchy
    _field_info_class = libytFieldInfo
    # _field_info_base_class = object
    _dataset_type = (
        "libyt"  # must set here since __init__() does not know dataset_type when calling it
    )
    _debug = False  # debug mode for libyt (not supported yet), some checks are in libyt C-library. Cannot open yet. Future use.
    libyt = None

    def __init__(self, units_override=None, unit_system="cgs"):
        # nothing to do if initialization has been done
        if self.libyt is not None:
            return

        # load the libyt module
        self.libyt = self._obtain_libyt()

        # get the target frontend (precisely speaking, the target Dataset subclass)
        # and set fluid type and fields accordingly
        self._code_frontend = self.libyt.param_yt["frontend"].lower()
        found_frontend = False
        for name in yt.frontends.api._frontends:
            if self._code_frontend == name.lower():
                # Import frontend dataset
                # The naming rule of dataset class is XXXDataset, XXX doesn't necessarily need to be all capital
                frontend = importlib.import_module(f"yt.frontends.{name.lower()}.api")
                frontend_dataset = None
                for _ in dir(frontend):
                    if _.endswith("Dataset"):
                        frontend_dataset = getattr(frontend, _)

                # Borrow frontend's field info
                self._field_info_class = frontend_dataset._field_info_class
                self.fluid_types += (self._code_frontend,)
                found_frontend = True

                # Read from param_yt['field_list'] and update the XXXDataset._field_info_class.
                # This action accumulates in each round, because we change the class
                # static variable in XXXDataset._field_info_class.
                if "field_list" in self.libyt.param_yt:
                    field_list = self.libyt.param_yt["field_list"]
                    known_other_fields = list(self._field_info_class.known_other_fields)
                    for field_name in field_list.keys():
                        # Step1 : Check if field_name is already inside known_other_fields
                        field_exist = False
                        for index in range(len(known_other_fields)):
                            if field_name == known_other_fields[index][0]:
                                field_exist = True
                                # Step2 : If field_name exists, append alias names one by one if it's not in the name
                                # list yet
                                for name_alias in field_list[field_name]["attribute"][1]:
                                    if name_alias not in known_other_fields[index][1][1]:
                                        known_other_fields[index][1][1].append(name_alias)
                                break
                        # Step2 : If field_name doesn't exist in known_other_fields, add a new field to add_fields list
                        if field_exist is False:
                            known_other_fields.append(
                                (field_name, tuple(field_list[field_name]["attribute"]))
                            )

                    # Step3 : convert it back to tuple
                    self._field_info_class.known_other_fields = tuple(known_other_fields)
                else:
                    mylog.debug("No self.libyt.param_yt['field_list'].")

                # Read from param_yt['particle_list'] and update XXXDataset._field_info_class.
                if "particle_list" in self.libyt.param_yt:
                    particle_list = self.libyt.param_yt["particle_list"]
                    known_particle_fields = list(self._field_info_class.known_particle_fields)
                    for ptype in particle_list.keys():
                        attribute = particle_list[ptype]["attribute"]
                        for particle in attribute.keys():
                            par_exist = False
                            for index in range(len(known_particle_fields)):
                                if particle == known_particle_fields[index][0]:
                                    par_exist = True
                                    for name_alias in attribute[particle][1]:
                                        if name_alias not in known_particle_fields[index][1][1]:
                                            known_particle_fields[index][1][1].append(name_alias)
                                    break

                            if par_exist is False:
                                known_particle_fields.append((particle, tuple(attribute[particle])))

                    self._field_info_class.known_particle_fields = tuple(known_particle_fields)
                else:
                    mylog.debug("No self.libyt.param_yt['particle_list'].")

                break
        if found_frontend is False:
            # We assume that user's code has corresponding yt frontend, if not, terminate yt.
            raise NotImplementedError(
                "libyt set frontend = %s, cannot find the code frontend [ %sDataset ] in yt."
                % (self.libyt.param_yt["frontend"], self.libyt.param_yt["frontend"].upper())
            )

        mylog.info("libyt: code dataset       = libytDataset")
        mylog.info(f"libyt: FieldInfo subclass = {self._field_info_class}")
        mylog.info(f"libyt: fluid type         = {self._code_frontend}")

        Dataset.__init__(
            self,
            self.libyt.param_yt["fig_basename"],
            self._dataset_type,
            units_override=units_override,
            unit_system=unit_system,
        )

    def _set_code_unit_attributes(self):
        # currently libyt assumes cgs
        setdefaultattr(self, "length_unit", self.quan(self.libyt.param_yt["length_unit"], "cm"))
        setdefaultattr(self, "mass_unit", self.quan(self.libyt.param_yt["mass_unit"], "g"))
        setdefaultattr(self, "time_unit", self.quan(self.libyt.param_yt["time_unit"], "s"))
        setdefaultattr(
            self, "magnetic_unit", self.quan(self.libyt.param_yt["magnetic_unit"], "gauss")
        )
        setdefaultattr(
            self, "velocity_unit", self.quan(self.libyt.param_yt["velocity_unit"], "cm/s")
        )

    def _parse_parameter_file(self):
        # dataset identifier
        self.unique_identifier = time.time()

        # user-specific parameters
        self.parameters.update(self.libyt.param_user)

        # yt-specific parameters
        param_yt = self.libyt.param_yt
        self.current_time = param_yt["current_time"]
        self.dimensionality = param_yt["dimensionality"]
        self.refine_by = param_yt["refine_by"]
        self.cosmological_simulation = param_yt["cosmological_simulation"]
        self.current_redshift = param_yt["current_redshift"]
        self.omega_matter = param_yt["omega_matter"]
        self.omega_lambda = param_yt["omega_lambda"]
        self.hubble_constant = param_yt["hubble_constant"]

        # vectors are stored as tuples in libyt and must be converted to NumPy arrays
        self.domain_left_edge = np.asarray(param_yt["domain_left_edge"])
        self.domain_right_edge = np.asarray(param_yt["domain_right_edge"])
        self.domain_dimensions = np.asarray(param_yt["domain_dimensions"])
        self._periodicity = (
            bool(param_yt["periodicity"][0] == 1),
            bool(param_yt["periodicity"][1] == 1),
            bool(param_yt["periodicity"][2] == 1),
        )

        # particle types
        if "particle_list" in param_yt:
            self.particle_types = tuple(param_yt["particle_list"].keys())
        else:
            self.particle_types = ()

        # Load code specific parameters
        for key in self.libyt.param_user.keys():
            if hasattr(self, key):
                mylog.info(
                    "Overwrite existing attribute self.%s = %s in class libytDataset",
                    key,
                    getattr(self, key),
                )

            setattr(self, key, self.libyt.param_user[key])
            mylog.info(
                "Set attribute self.%s = %s in class libytDataset.",
                key,
                self.libyt.param_user[key],
            )

    @staticmethod
    def _obtain_libyt():
        import libyt

        from ._version import __version__

        libyt_version = f"{libyt.libyt_info['version'][0]}.{libyt.libyt_info['version'][1]}.{libyt.libyt_info['version'][2]}"
        mylog.info(f"libyt    version = {libyt_version}")
        mylog.info(f"yt_libyt version = {__version__}")

        if libyt.libyt_info["version"][0] > 0:
            mylog.error(
                "libyt version %s is not compatible with yt_libyt version %s. "
                "Expecting libyt >=0.1.0,<1.0.0." % (libyt_version, __version__)
            )
            raise ValueError(
                "libyt version %s is not compatible with yt_libyt version %s. "
                "Expecting libyt >=0.1.0,<1.0.0." % (libyt_version, __version__)
            )

        return libyt

    @classmethod
    def _is_valid(cls, *args, **kwargs):
        # always return false since libyt is used for inline analysis only
        return False
