"""
libyt-specific IO functions



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.utilities.io_handler import BaseIOHandler
from yt.funcs import mylog
from yt.geometry.selection_routines import AlwaysSelector


class IOHandlerlibyt(BaseIOHandler):
    _particle_reader = False
    _dataset_type    = "libyt"

    def __init__(self, ds):
        super(IOHandlerlibyt, self).__init__(ds)
        import libyt
        self.libyt        = libyt
        self.ds           = ds
        self.grid_data    = libyt.grid_data
        self.param_yt     = libyt.param_yt
        self.hierarchy    = libyt.hierarchy
        self._field_dtype = "float64"

###     ghost_zones != 0 is not supported yet
#       self.my_slice = (slice(ghost_zones,-ghost_zones),
#                        slice(ghost_zones,-ghost_zones),
#                        slice(ghost_zones,-ghost_zones))

    def _read_particle_coords(self, chunks, ptf):
        chunks = list(chunks)
        # TODO: Support get non-local grid
        #  (1) Collect all the particle attribute in a dictionary, since we need to make
        #      all the rank in the same RMA epoch.

        for chunk in chunks:
            for g in chunk.objs:
                # if grid_particle_count, which is sum of all particle number
                # in that grid is zero, continue
                if self.hierarchy['grid_particle_count'][g.id] == 0:
                    continue

                # else, fetch the position x/y/z of particle by ptype
                for ptype in ptf.keys():
                    coor_label = self.param_yt['particle_list'][ptype]['particle_coor_label']
                    if None in coor_label:
                        raise ValueError("coor_x, coor_y, coor_z label not set!")
                    x = self.libyt.get_attr(g.id, ptype, coor_label[0])
                    y = self.libyt.get_attr(g.id, ptype, coor_label[1])
                    z = self.libyt.get_attr(g.id, ptype, coor_label[2])

                    # g.id ptype particle number is 0, libyt.get_attr will return None
                    if x is None or y is None or z is None:
                        continue
                    else:
                        yield ptype, (x, y, z)

    def _read_particle_fields(self, chunks, ptf, selector):
        chunks = list(chunks)
        # TODO: Support get non-local grid
        #  (1) Collect all the particle attribute in a dictionary, since we need to make
        #      all the rank in the same RMA epoch.
        for chunk in chunks:
            for g in chunk.objs:
                # if grid_particle_count, which is sum of all particle number
                # in that grid is zero, continue
                if self.hierarchy['grid_particle_count'][g.id] == 0:
                    continue

                # else, fetch the position x/y/z of particle by ptype
                for ptype in ptf.keys():
                    coor_label = self.param_yt['particle_list'][ptype]['particle_coor_label']
                    if None in coor_label:
                        raise ValueError("coor_x, coor_y, coor_z label not set!")
                    x = self.libyt.get_attr(g.id, ptype, coor_label[0])
                    y = self.libyt.get_attr(g.id, ptype, coor_label[1])
                    z = self.libyt.get_attr(g.id, ptype, coor_label[2])

                    # g.id ptype particle number is 0, libyt.get_attr will return None
                    if x is None or y is None or z is None:
                        continue

                    mask = selector.select_points(x, y, z, 0.0)
                    if mask is None:
                        continue

                    for field in ptf[ptype]:
                        data = self.libyt.get_attr(g.id, ptype, field)
                        # if ptype particle num in grid g.id = 0, get_attr will return None.
                        # then we shall continue the loop
                        if data is None:
                            continue
                        else:
                            yield (ptype, field), data[mask]

    def _read_chunk_data(self, chunk, fields):
        # TODO: The suite hasn't been tested yet.
        #       Although it's be use for caching, I wonder do libyt need this.
        #       Since we don't need to load data from file. Although we do need
        #       to get data from remote rank.
        rv = {}
        if len(chunk.objs) == 0:
            return rv
        for g in chunk.objs:
            rv[g.id] = {}

        # Split into particles and non-particles
        fluid_fields, particle_fields = [], []
        for ftype, fname in fields:
            if ftype in self.ds.particle_types:
                particle_fields.append((ftype, fname))
            else:
                fluid_fields.append((ftype, fname))

        # Read particle data
        if len(particle_fields) > 0:
            selector = AlwaysSelector(self.ds)
            rv.update(self._read_particle_selection([chunk], selector, particle_fields))

        # If no more fluid fields to read, return rv. Or else, read fluid fields
        if len(fluid_fields) == 0:
            return rv

        for field in fluid_fields:
            ftype, fname = field
            for g in chunk.objs:
                rv[g.id][field] = self._get_field_from_libyt(g, fname)
        return rv

    def _read_fluid_selection(self, chunks, selector, fields, size):
        rv = {}
        chunks = list(chunks)

        if selector.__class__.__name__ == "GridSelector":
            if not (len(chunks) == len(chunks[0].objs) == 1):
                raise RuntimeError("class IOHandlerlibyt, def _read_fluid_selection, selector == GridSelector, "
                                   "chunk to be read not equal to 1.")
            g = chunks[0].objs[0]
            for ftype, fname in fields:
                rv[(ftype, fname)] = self._get_field_from_libyt(g, fname)
            return rv

        if size is None:
            size = sum((g.count(selector) for chunk in chunks for g in chunk.objs))

        for field in fields:
            rv[field] = np.empty(size, dtype=self._field_dtype)

        ng = sum(len(c.objs) for c in chunks)
        mylog.debug("Reading %s cells of %s fields in %s grids",
                    size, [f2 for f1, f2 in fields], ng)

        for ftype, fname in fields:
            mylog.debug("ftype = %s" % ftype)
            mylog.debug("fname = %s" % fname)

        # TODO: Support get non-local grid,
        #  (1) Maybe I should make _get_data_from_libyt method return in a group of needed grid, not just a single grid.
        #  (2) Inside _get_data_from_libyt, there should be two methods one returns local grids, the other returns
        #      non-local grids.

        # Distinguish local and non-local grid
        local_id, to_prepare, nonlocal_id, nonlocal_rank = self._distinguish_nonlocal_grids(chunks)

        # Get local grid
        for field in fields:
            offset = 0
            ftype, fname = field
            for chunk in chunks:
                for g in chunk.objs:
                    data_view = self._get_field_from_libyt(g, fname)
                    offset += g.select(selector, data_view, rv[field], offset)
            assert (offset == size)
        return rv

    def _distinguish_nonlocal_grids(self, chunks):
        # Split local and non-local grids.
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        myrank = comm.Get_rank()

        local_id = []
        nonlocal_id = []
        nonlocal_rank = []

        for chunk in chunks:
            for g in chunk.objs:
                if g.MPI_rank != myrank:
                    nonlocal_id.append(g.id)
                    nonlocal_rank.append(g.MPI_rank)
                else:
                    local_id.append(g.id)
        num_nonlocal_grids = len(nonlocal_id)

        mylog.debug("local grid = %s" % local_id)
        mylog.debug("num nonlocal grid  = %s" % num_nonlocal_grids)
        mylog.debug("nonlocal grid id   = %s" % nonlocal_id)
        mylog.debug("nonlocal grid rank = %s" % nonlocal_rank)

        # Gather all non-local grids in each rank.
        sendcounts = comm.gather(num_nonlocal_grids, root=0)
        sendcounts = comm.bcast(sendcounts, root=0)
        sendbuf = np.asarray(nonlocal_id)
        recvbuf = np.empty(sum(sendcounts), dtype=int)
        comm.Gatherv(sendbuf=sendbuf, recvbuf=(recvbuf, sendcounts), root=0)
        comm.Bcast(recvbuf, root=0)

        mylog.debug("sendcounts = %s" % sendcounts)
        mylog.debug("all non-local grids = %s" % recvbuf)

        # Get grid id that this rank has to prepare.
        proc_num = self.hierarchy["proc_num"][:, 0]
        index = np.argwhere(proc_num[recvbuf] == myrank)
        to_prepare = list(np.unique(recvbuf[index]))

        mylog.debug("to_prepare = %s" % to_prepare)

        return local_id, to_prepare, nonlocal_id, nonlocal_rank

    def _get_remote_field_from_libyt(self, grids, fields):
        # Wrapper for the RMA operation at libyt C library code.
        pass

    def _get_remote_particle_from_libyt(self, grids, ptf):
        # Counter-part of _get_remote_field_from_libyt for supporting particles.
        pass

    def _get_field_from_libyt(self, grid, fname):
        # This method is to get the local grid data.
        field_list = self.param_yt["field_list"]
        if field_list[fname]["field_define_type"] == "cell-centered":
            data_convert = self.grid_data[grid.id][fname]
            assert data_convert is not None, "This MPI rank does not have grid id [%s]." % grid.id
        elif field_list[fname]["field_define_type"] == "face-centered":
            # convert to cell-centered
            data_temp = self.grid_data[grid.id][fname]
            assert data_temp is not None, "This MPI rank does not have grid id [%s]." % grid.id
            grid_dim = self.hierarchy["grid_dimensions"][grid.id]
            if field_list[fname]["swap_axes"] is True:
                grid_dim = np.flip(grid_dim)
            axis = np.argwhere(grid_dim != data_temp.shape)
            assert len(axis) == 1, \
                "Field [ %s ] is not a face-centered data, " \
                "grid_dimensions = %s, field data dimensions = %s, considering swap_axes" % (fname, grid_dim, (data_temp.shape,))
            assert data_temp.shape[axis[0, 0]] - 1 == grid_dim[axis[0, 0]], \
                "Field [ %s ] is not a face-centered data, " \
                "grid_dimensions = %s, field data dimensions = %s, considering swap_axes" % (fname, grid_dim, (data_temp.shape,))
            if axis == 0:
                data_convert = 0.5 * (data_temp[:-1, :, :] + data_temp[1:, :, :])
            elif axis == 1:
                data_convert = 0.5 * (data_temp[:, :-1, :] + data_temp[:, 1:, :])
            elif axis == 2:
                data_convert = 0.5 * (data_temp[:, :, :-1] + data_temp[:, :, 1:])
        elif field_list[fname]["field_define_type"] == "derived_func":
            data_convert = self.libyt.derived_func(grid.id, fname)
        else:
            # Since we only supports "cell-centered", "face-centered", "derived_func" tags for now
            # Raise an error if enter this block.
            raise ValueError("libyt does not have field_define_type [ %s ]" %
                             (field_list[fname]["field_define_type"]))

        # Swap axes or not, then return
        if field_list[fname]["swap_axes"] is True:
            return data_convert.swapaxes(0, 2)
        else:
            return data_convert
