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

import os
import stat
import time
import numpy as np
import weakref

from yt.funcs import \
    mylog, \
    setdefaultattr
from yt.data_objects.index_subobjects.grid_patch import \
    AMRGridPatch
from yt.geometry.grid_geometry_handler import \
    GridIndex
from yt.data_objects.static_output import \
    Dataset
from yt.utilities.file_handler import \
    HDF5FileHandler
from yt.utilities.parameter_file_storage import \
    output_type_registry
from .fields import libytFieldInfo


class libytGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, index, level):
        AMRGridPatch.__init__(self, id,
                              filename=None,
                              index=index)
        self.Parent = None
        self.Children = []
        self.Level = level

    def __repr__(self):
        return 'libytGrid_%09i (dimension = %s)' % (self.id, self.ActiveDimensions)


class libytHierarchy(GridIndex):
    grid = libytGrid
    libyt = None

    ### NOT SURE ABOUT THIS OPTION
    #   _preload_implemented = True

    def __init__(self, ds, dataset_type='libyt'):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        self.directory = os.getcwd()
        self.float_type = 'float64'
        self.libyt = self.dataset.libyt

        GridIndex.__init__(self, ds, dataset_type)

    def _detect_output_fields(self):
        # assuming all grids have the same fields
        gid = 0
        self.field_list = [(self.dataset._fluid_type, v)
                           for v in self.libyt.grid_data[gid].keys()]

        # assuming all particles have the same fields

    #       if self._group_particle is not None:
    #           self.field_list += [ ('io', v) for v in self._group_particle.keys() ]

    def _count_grids(self):
        self.num_grids = self.libyt.param_yt['num_grids']

    def _parse_index(self):
        # hierarchy information of all grids
        hierarchy = self.libyt.hierarchy

        # **copy** data from hierarchy since they may have different data types than self.grid_xxx
        # also note that all 1D arrays below are actually declared as [:][1]
        self.grid_left_edge[:] = hierarchy['grid_left_edge'][:]
        self.grid_right_edge[:] = hierarchy['grid_right_edge'][:]
        self.grid_dimensions[:] = hierarchy['grid_dimensions'][:]
        self.grid_levels[:] = hierarchy['grid_levels'][:]
        self.grid_particle_count[:] = hierarchy['grid_particle_count'][:]
        self.max_level = self.grid_levels.max()

        # allocate all grid objects
        self.grids = np.empty(self.num_grids, dtype='object')
        for i in range(self.num_grids):
            self.grids[i] = self.grid(i, self, self.grid_levels.flat[i])

        # calculate the starting particle indices for each grid (starting from 0)
        # --> note that the last element must store the total number of particles
        #    (see _read_particle_coords and _read_particle_fields in io.py)
        # self._particle_indices = np.zeros(self.num_grids + 1, dtype='int64')
        # np.add.accumulate(self.grid_particle_count.squeeze(), out=self._particle_indices[1:])

    def _populate_grid_objects(self):
        # must flat it since 'grid_parent_id' has the dimension [num_grids][1]
        parent_list = self.libyt.hierarchy['grid_parent_id'].flat

        for gid in range(self.num_grids):
            grid = self.grids[gid]
            parent_id = parent_list[gid]

            # set up the parent-children relationship
            if parent_id >= 0:
                parent_grid = self.grids[parent_id]
                grid.Parent = parent_grid
                parent_grid.Children.append(grid)

            # set up other grid attributes
            grid._prepare_grid()
            grid._setup_dx()

        # validate the parent-children relationship in the debug mode
        if self.dataset._debug:
            self._validate_parent_children_relasionship()

    #   # for _debug mode only (not implemented yet)
    def _validate_parent_children_relasionship(self):
        pass


#       mylog.info('Validating the parent-children relationship ...')
#       mylog.info('Check passed')

### NOT SURE IF WE NEED THIS (see EnzoHierarchyInMemory class)
#   def _chunk_io(self, dobj, cache = True, local_only = False):
#       gfiles = defaultdict(list)
#       gobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
#       for g in gobjs:
#           gfiles[g.filename].append(g)
#       for fn in sorted(gfiles):
#           if local_only:
#               gobjs = [g for g in gfiles[fn] if g.proc_num == self.comm.rank]
#               gfiles[fn] = gobjs
#           gs = gfiles[fn]
#           count = self._count_selection(dobj, gs)
#           yield YTDataChunk(dobj, "io", gs, count, cache = cache)


class libytDataset(Dataset):
    _index_class = libytHierarchy
    _field_info_class = libytFieldInfo
    _dataset_type = 'libyt'  # must set here since __init__() does not know dataset_type when calling it
    _debug = False  # debug mode for libyt (not supported yet)
    libyt = None

    def __init__(self,
                 units_override=None,
                 unit_system="cgs"):

        # nothing to do if initialization has been done
        if self.libyt is not None: return

        # load the libyt module
        self.libyt = self._obtain_libyt()

        # get the target frontend (precisely speaking, the target Dataset subclass)
        # and set fluid type and fields accordingly
        self._code_frontend = self.libyt.param_yt['frontend'].lower()
        for name, cls in output_type_registry.items():

            mylog.debug("#FLAG#")
            mylog.debug("yt/frontends/libyt/data_structures.py (class libytDataset, def __init__)")
            mylog.debug("name = %s", name)
            mylog.debug("cls = %s", cls)

            if self._code_frontend == name[0:-len('Dataset')].lower():
                mylog.debug("self.fluid_types = %s", self.fluid_types)

                self._code_dataset = name
                self._field_info_class = cls._field_info_class

                # No _fluid_type in yt-3.6.0 at all
                # self._fluid_type       = cls._fluid_type
                # self.fluid_types      += ( self._fluid_type, )

                self._fluid_type = self._code_frontend
                self.fluid_types += (self._fluid_type,)

                break
        else:
            mylog.warning('Cannot find the code frontend [%s]' % self.libyt.param_yt['frontend'])
            # set various attributes to libyt itself
            self._code_dataset = self.__class__.__name__
            self._field_info_class = libytFieldInfo
            self._fluid_type = 'libyt'
            self.fluid_types += (self._fluid_type,)

        mylog.info('libyt: code dataset       = %s' % self._code_dataset)
        mylog.info('libyt: FieldInfo subclass = %s' % self._field_info_class)
        mylog.info('libyt: fluid type         = %s' % self._fluid_type)

        Dataset.__init__(self, "libytHasNoParameterFile", self._dataset_type,
                         units_override=units_override,
                         unit_system=unit_system)

    def _set_code_unit_attributes(self):
        # currently libyt assumes cgs
        setdefaultattr(self, 'length_unit', self.quan(self.libyt.param_yt['length_unit'], 'cm'))
        setdefaultattr(self, 'mass_unit', self.quan(self.libyt.param_yt['mass_unit'], 'g'))
        setdefaultattr(self, 'time_unit', self.quan(self.libyt.param_yt['time_unit'], 's'))

    def _parse_parameter_file(self):
        # dataset identifier
        self.unique_identifier = time.time()

        # user-specific parameters
        self.parameters.update(self.libyt.param_user)

        # yt-specific parameters
        param_yt = self.libyt.param_yt
        self.basename = param_yt['fig_basename']
        self.current_time = param_yt['current_time']
        self.dimensionality = param_yt['dimensionality']
        self.refine_by = param_yt['refine_by']
        self.cosmological_simulation = param_yt['cosmological_simulation']
        self.current_redshift = param_yt['current_redshift']
        self.omega_matter = param_yt['omega_matter']
        self.omega_lambda = param_yt['omega_lambda']
        self.hubble_constant = param_yt['hubble_constant']

        # vectors are stored as tuples in libyt and must be converted to NumPy arrays
        self.domain_left_edge = np.asarray(param_yt['domain_left_edge'])
        self.domain_right_edge = np.asarray(param_yt['domain_right_edge'])
        self.domain_dimensions = np.asarray(param_yt['domain_dimensions'])
        self.periodicity = np.asarray(param_yt['periodicity'])

        # Just to make example/example run
        if (self._code_frontend == "gamer"):
            self.mhd = 0

    def _obtain_libyt(self):
        import libyt
        return libyt

    @classmethod
    def _is_valid(cls, *args, **kwargs):
        # always return false since libyt is used for inline analysis only
        return False

