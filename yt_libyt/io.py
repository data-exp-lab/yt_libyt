"""
libyt-specific IO functions



"""

# -----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

import numpy as np
from yt.funcs import mylog
from yt.geometry.selection_routines import AlwaysSelector
from yt.utilities.io_handler import BaseIOHandler


class libytIOHandler(BaseIOHandler):
    _particle_reader = False
    _dataset_type = "libyt"

    def __init__(self, ds):
        super(libytIOHandler, self).__init__(ds)
        import libyt

        self.libyt = libyt
        self.ds = ds
        self.grid_data = libyt.grid_data
        self.param_yt = libyt.param_yt
        self.hierarchy = libyt.hierarchy
        self._field_dtype = "float64"
        self.myrank = libytIOHandler._get_my_rank()

    def _read_particle_coords(self, chunks, ptf):
        chunks = list(chunks)

        # Get position (coordinate) label.
        ptf_new = {}
        for ptype in ptf.keys():
            coor_label = self.param_yt["particle_list"][ptype]["particle_coor_label"]
            if None in coor_label:
                raise ValueError("Particle label representing postion X/Y/Z not set!")
            ptf_new[ptype] = list(coor_label)

        # Get remote data.
        nonlocal_data = self._prepare_remote_particle_from_libyt(chunks, ptf_new)

        # Get index offset and particle data dict
        index_offset = self.param_yt["index_offset"]
        particle_data = self.libyt.particle_data

        for chunk in chunks:
            for g in chunk.objs:
                # fetch the position x/y/z of particle by ptype
                for ptype in ptf.keys():
                    # Get particle count in ptype, continue if it is zero
                    index_label = self.param_yt["particle_list"][ptype]["label"]
                    if self.hierarchy["par_count_list"][g.id - index_offset][index_label] == 0:
                        continue

                    coor_label = self.param_yt["particle_list"][ptype]["particle_coor_label"]
                    if g.MPI_rank == self.myrank:
                        if (
                            g.id in particle_data
                            and ptype in particle_data[g.id]
                            and coor_label[0] in particle_data[g.id][ptype]
                        ):
                            x = particle_data[g.id][ptype][coor_label[0]]
                        else:
                            x = self.libyt.get_particle(g.id, ptype, coor_label[0])

                        if (
                            g.id in particle_data
                            and ptype in particle_data[g.id]
                            and coor_label[1] in particle_data[g.id][ptype]
                        ):
                            y = particle_data[g.id][ptype][coor_label[1]]
                        else:
                            y = self.libyt.get_particle(g.id, ptype, coor_label[1])

                        if (
                            g.id in particle_data
                            and ptype in particle_data[g.id]
                            and coor_label[2] in particle_data[g.id][ptype]
                        ):
                            z = particle_data[g.id][ptype][coor_label[2]]
                        else:
                            z = self.libyt.get_particle(g.id, ptype, coor_label[2])
                    else:
                        x = nonlocal_data[g.id][ptype][coor_label[0]]
                        y = nonlocal_data[g.id][ptype][coor_label[1]]
                        z = nonlocal_data[g.id][ptype][coor_label[2]]

                    # g.id ptype particle number is 0, libyt.get_particle will return None,
                    # It will not happen unless something went wrong when passing particle count to libyt.
                    if x is None or y is None or z is None:
                        raise ValueError("Particle position should not be None.")
                    else:
                        yield ptype, (x, y, z)

    def _read_particle_fields(self, chunks, ptf, selector):
        chunks = list(chunks)

        # Get position (coordinate) label and append particle attribute to get after them.
        ptf_new = {}
        for ptype in ptf.keys():
            coor_label = self.param_yt["particle_list"][ptype]["particle_coor_label"]
            if None in coor_label:
                raise ValueError("Particle label representing postion X/Y/Z not set!")
            ptf_new[ptype] = list(coor_label)
            for field in ptf[ptype]:
                ptf_new[ptype].append(field)
            ptf_new[ptype] = set(ptf_new[ptype])

        # Get remote data.
        nonlocal_data = self._prepare_remote_particle_from_libyt(chunks, ptf_new)

        # Get index offset and particle data
        index_offset = self.param_yt["index_offset"]
        particle_data = self.libyt.particle_data

        for chunk in chunks:
            for g in chunk.objs:
                # fetch particle data.
                for ptype in ptf.keys():
                    # get particle count in ptype, continue if it is zero
                    index_label = self.param_yt["particle_list"][ptype]["label"]
                    if self.hierarchy["par_count_list"][g.id - index_offset][index_label] == 0:
                        continue

                    # fetch the position x/y/z of particle by ptype
                    coor_label = self.param_yt["particle_list"][ptype]["particle_coor_label"]
                    if None in coor_label:
                        raise ValueError("Particle label representing postion X/Y/Z not set!")
                    if g.MPI_rank == self.myrank:
                        if (
                            g.id in particle_data
                            and ptype in particle_data[g.id]
                            and coor_label[0] in particle_data[g.id][ptype]
                        ):
                            x = particle_data[g.id][ptype][coor_label[0]]
                        else:
                            x = self.libyt.get_particle(g.id, ptype, coor_label[0])

                        if (
                            g.id in particle_data
                            and ptype in particle_data[g.id]
                            and coor_label[1] in particle_data[g.id][ptype]
                        ):
                            y = particle_data[g.id][ptype][coor_label[1]]
                        else:
                            y = self.libyt.get_particle(g.id, ptype, coor_label[1])

                        if (
                            g.id in particle_data
                            and ptype in particle_data[g.id]
                            and coor_label[2] in particle_data[g.id][ptype]
                        ):
                            z = particle_data[g.id][ptype][coor_label[2]]
                        else:
                            z = self.libyt.get_particle(g.id, ptype, coor_label[2])
                    else:
                        x = nonlocal_data[g.id][ptype][coor_label[0]]
                        y = nonlocal_data[g.id][ptype][coor_label[1]]
                        z = nonlocal_data[g.id][ptype][coor_label[2]]

                    # g.id ptype particle number is 0, libyt.get_particle will return None.
                    # It will not happen unless something went wrong when passing particle count to libyt.
                    if x is None or y is None or z is None:
                        raise ValueError("Particle position should not be None.")

                    mask = selector.select_points(x, y, z, 0.0)
                    if mask is None:
                        continue

                    for field in ptf[ptype]:
                        if g.MPI_rank == self.myrank:
                            if (
                                g.id in particle_data
                                and ptype in particle_data[g.id]
                                and field in particle_data[g.id][ptype]
                            ):
                                data = particle_data[g.id][ptype][field]
                            else:
                                data = self.libyt.get_particle(g.id, ptype, field)
                        else:
                            data = nonlocal_data[g.id][ptype][field]

                        # if ptype particle num in grid g.id = 0, get_particle will return None.
                        # It will not happen unless something went wrong when passing particle count to libyt.
                        if data is None:
                            raise ValueError("Particle data should not be None.")
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

        # Prepare nonlocal data
        nonlocal_data = self._prepare_remote_field_from_libyt([chunk], fields)

        # TODO: Bug, rv should allocate a new buffer.
        for field in fluid_fields:
            ftype, fname = field
            for g in chunk.objs:
                if g.MPI_rank == self.myrank:
                    rv[g.id][field] = self._get_field_from_libyt(g, fname)
                else:
                    rv[g.id][field] = self._get_field_from_libyt(
                        g, fname, nonlocal_data=nonlocal_data
                    )
        return rv

    def _read_fluid_selection(self, chunks, selector, fields, size):
        rv = {}
        chunks = list(chunks)

        # Prepare nonlocal data
        nonlocal_data = self._prepare_remote_field_from_libyt(chunks, fields)

        # TODO: Allocate buffer for rv, don't make rv point directly to simulation data buffer.
        # if selector.__class__.__name__ == "GridSelector":
        #     if not (len(chunks) == len(chunks[0].objs) == 1):
        #         raise RuntimeError("class libytIOHandler, def _read_fluid_selection, selector == GridSelector, "
        #                            "chunk to be read not equal to 1.")
        #     g = chunks[0].objs[0]
        #     for ftype, fname in fields:
        #         if g.MPI_rank == self.myrank:
        #             rv[(ftype, fname)] = self._get_field_from_libyt(g, fname)
        #         else:
        #             rv[(ftype, fname)] = self._get_field_from_libyt(g, fname, nonlocal_data=nonlocal_data)
        #     return rv

        if size is None:
            size = sum((g.count(selector) for chunk in chunks for g in chunk.objs))

        for field in fields:
            rv[field] = np.empty(size, dtype=self._field_dtype)

        ng = sum(len(c.objs) for c in chunks)
        mylog.debug(
            "Reading %s cells of %s fields in %s grids", size, [f2 for f1, f2 in fields], ng
        )

        # Get grid data
        for field in fields:
            offset = 0
            ftype, fname = field
            for chunk in chunks:
                for g in chunk.objs:
                    if g.MPI_rank == self.myrank:
                        data_view = self._get_field_from_libyt(g, fname)
                    else:
                        data_view = self._get_field_from_libyt(
                            g, fname, nonlocal_data=nonlocal_data
                        )
                    offset += g.select(selector, data_view, rv[field], offset)
            assert offset == size
        return rv

    @staticmethod
    def _get_my_rank():
        import libyt

        if libyt.libyt_info["SERIAL_MODE"] is False:
            try:
                from mpi4py import MPI

                comm = MPI.COMM_WORLD
                return comm.Get_rank()
            except ImportError:
                raise ImportError("Need mpi4py in parallel mode (SERIAL_MODE = false)")
        else:
            return 0

    def _distinguish_nonlocal_grids(self, chunks):

        if self.libyt.libyt_info["SERIAL_MODE"] is True:
            # we don't need rma
            return False, [], [], []

        # Split local and non-local grids.
        from mpi4py import MPI

        comm = MPI.COMM_WORLD

        local_id = []
        nonlocal_id = []
        nonlocal_rank = []

        for chunk in chunks:
            for g in chunk.objs:
                if g.MPI_rank != self.myrank:
                    nonlocal_id.append(g.id)
                    nonlocal_rank.append(g.MPI_rank)
                else:
                    local_id.append(g.id)
        num_nonlocal_grids = len(nonlocal_id)

        # Gather all non-local grids in each rank.
        sendcounts = comm.gather(num_nonlocal_grids, root=0)
        sendcounts = comm.bcast(sendcounts, root=0)
        sendbuf = np.asarray(nonlocal_id)
        recvbuf = np.empty(sum(sendcounts), dtype=int)
        comm.Gatherv(sendbuf=sendbuf, recvbuf=(recvbuf, sendcounts), root=0)
        comm.Bcast(recvbuf, root=0)

        # Get grid id that this rank has to prepare, grid id doesn't have to be 0-indexed
        index_offset = self.param_yt["index_offset"]
        proc_num = self.hierarchy["proc_num"][:, 0]
        index = np.argwhere(proc_num[recvbuf - index_offset] == self.myrank)
        to_prepare = list(np.unique(recvbuf[index]))

        # Determine whether we should call for libyt C extend method for RMA operation.
        if sum(sendcounts) != 0:
            rma = True
        else:
            rma = False

        return rma, to_prepare, nonlocal_id, nonlocal_rank

    def _prepare_remote_field_from_libyt(self, chunks, fields):
        # Wrapper for the RMA operation at libyt C library code.
        # Each rank must call this method, in order to get nonlocal grids.

        # Distinguish local and non-local grid, and what should this rank prepared.
        rma, to_prepare, nonlocal_id, nonlocal_rank = self._distinguish_nonlocal_grids(chunks)

        if rma is True:
            # Encode field name to UTF-8
            fname_list = []
            for _ftype, fname in fields:
                fname_list.append(fname.encode(encoding="UTF-8", errors="strict"))
            fname_list = sorted(set(fname_list))

            # Get nonlocal_data, libyt will perform RMA operation in this step.
            # Every rank must call this libyt method.
            mylog.debug("Getting nonlocal data through libyt ...")
            nonlocal_data = self.libyt.get_field_remote(
                fname_list,
                len(fname_list),
                to_prepare,
                len(to_prepare),
                nonlocal_id,
                nonlocal_rank,
                len(nonlocal_id),
            )
        else:
            nonlocal_data = None

        return nonlocal_data

    def _prepare_remote_particle_from_libyt(self, chunks, ptf):
        # Wrapper for the RMA operation at libyt C library code.
        # For supporting particles. Each rank must call this method.

        # Distinguish local and non-local grid, and what should this rank prepared.
        rma, to_prepare, nonlocal_id, nonlocal_rank = self._distinguish_nonlocal_grids(chunks)

        if rma is True:
            # Filter out those who really has particles in their grid. Grid id doesn't have to be 0-indexed.
            # Since we aren't sure how many particle types will yt access, we check total particle counts.
            index_offset = self.param_yt["index_offset"]
            par_count = self.ds.index.grid_particle_count[:, 0]

            # If to_prepare is empty list, we don't need to filter out grids without particles
            if len(to_prepare) != 0:
                to_prepare = np.asarray(to_prepare)
                index = np.argwhere(par_count[to_prepare - index_offset] > 0)
                to_prepare = list(to_prepare[index].flatten())

            # If nonlocal_id / nonlocal_rank is empty list, we don't need to filter out grids without particles
            if len(nonlocal_id) != 0:
                nonlocal_id = np.asarray(nonlocal_id)
                index = np.argwhere(par_count[nonlocal_id - index_offset] > 0)
                nonlocal_id = list(nonlocal_id[index].flatten())
                nonlocal_rank = np.asarray(nonlocal_rank)
                nonlocal_rank = list(nonlocal_rank[index].flatten())

            # String inside ptf should be encoded in UTF-8, and attributes should be in list obj.
            ptf_c = {}
            for key in ptf.keys():
                attr_list = []
                for attr in ptf[key]:
                    attr_list.append(attr.encode(encoding="UTF-8", errors="strict"))
                ptype = key.encode(encoding="UTF-8", errors="strict")
                ptf_c[ptype] = sorted(set(attr_list))

            # Call libyt RMA
            mylog.debug("Getting nonlocal data through libyt ...")
            nonlocal_data = self.libyt.get_particle_remote(
                ptf_c,
                ptf_c.keys(),
                to_prepare,
                len(to_prepare),
                nonlocal_id,
                nonlocal_rank,
                len(nonlocal_id),
            )
        else:
            nonlocal_data = None

        return nonlocal_data

    def _remove_ghost_cells(self, ghost_cell, data):
        dim = data.ndim
        shape = data.shape
        if dim == 3:
            return data[
                ghost_cell[0] : (shape[0] - ghost_cell[1]),
                ghost_cell[2] : (shape[1] - ghost_cell[3]),
                ghost_cell[4] : (shape[2] - ghost_cell[5]),
            ]
        elif dim == 2:
            return data[
                ghost_cell[0] : (shape[0] - ghost_cell[1]),
                ghost_cell[2] : (shape[1] - ghost_cell[3]),
            ]
        elif dim == 1:
            return data[ghost_cell[0] : (shape[0] - ghost_cell[1])]
        else:
            mylog.error(
                "libyt does not support removing ghost cells for data with dimension %d." % dim
            )
            raise ValueError("libyt does not support removing ghost cells for this data.")

    def _convert_face_centered_to_cell_centered(self, fname, contiguous_in_x, grid_id, data):
        grid_dim = self.hierarchy["grid_dimensions"][grid_id - self.param_yt["index_offset"]][
            : self.param_yt["dimensionality"]
        ].copy()
        if contiguous_in_x is True:
            grid_dim = np.flip(grid_dim)
        axis = np.argwhere(grid_dim != data.shape).flatten()
        if len(axis) != 1 or data.shape[axis[0]] - 1 != grid_dim[axis[0]]:
            mylog.error(
                "Field [%s] in grid [%d] is not a face-centered data. "
                "Data grid has dim %s, but it cannot be converted to dim = %s"
                % (fname, grid_id, (data.shape,), grid_dim)
            )
            raise ValueError("Face-centered data dimension not match.")

        if self.param_yt["dimensionality"] == 3:
            if axis == 0:
                return 0.5 * (data[:-1, :, :] + data[1:, :, :])
            elif axis == 1:
                return 0.5 * (data[:, :-1, :] + data[:, 1:, :])
            else:
                return 0.5 * (data[:, :, :-1] + data[:, :, 1:])
        elif self.param_yt["dimensionality"] == 2:
            if axis == 0:
                return 0.5 * (data[:-1, :] + data[1:, :])
            else:
                return 0.5 * (data[:, :-1] + data[:, 1:])
        elif self.param_yt["dimensionality"] == 1:
            return 0.5 * (data[:-1] + data[1:])
        else:
            raise ValueError(
                "Unable to convert face-centered data to cell-centered data for dimensionality %d."
                % self.param_yt["dimensionality"]
            )

    def _get_field_from_libyt(self, grid, fname, nonlocal_data=None):
        # This method is to get the grid data.
        # If nonlocal_data is none, which means to get a local grid.
        # Otherwise, read the nonlocal data in nonlocal_data.
        field_list = self.param_yt["field_list"]
        ghost_cell = field_list[fname]["ghost_cell"]
        if field_list[fname]["field_type"] == "cell-centered":
            # Read data from grid_data, or nonlocal_data.
            # We don't create key-value pair if no data pass in from user.
            try:
                if nonlocal_data is None:
                    data_convert = self.grid_data[grid.id][fname]
                else:
                    data_convert = nonlocal_data[grid.id][fname]
            except Exception as err:
                mylog.error("%s: %s", type(err).__name__, str(err))
                mylog.error(
                    "Cannot get cell-centered grid [%s] data on MPI rank [%d]."
                    % (grid.id, grid.MPI_rank)
                )
                raise RuntimeError("libyt didn't get the data successfully.")

            # Remove ghost cell, and get my slice
            data_convert = self._remove_ghost_cells(ghost_cell, data_convert)

        elif field_list[fname]["field_type"] == "face-centered":
            # Read data from grid_data, or nonlocal_data.
            # We don't create key-value pair if no data pass in from user.
            try:
                if nonlocal_data is None:
                    data_temp = self.grid_data[grid.id][fname]
                else:
                    data_temp = nonlocal_data[grid.id][fname]
            except Exception as err:
                mylog.error("%s: %s", type(err).__name__, str(err))
                mylog.error(
                    "Cannot get face-centered grid [%s] data on MPI rank [%d]."
                    % (grid.id, grid.MPI_rank)
                )
                raise RuntimeError("libyt didn't get the data successfully.")

            # Remove ghost cell, and get my slice
            data_temp = self._remove_ghost_cells(ghost_cell, data_temp)

            # Convert to cell-centered
            data_convert = self._convert_face_centered_to_cell_centered(
                fname, field_list[fname]["contiguous_in_x"], grid.id, data_temp
            )

        elif field_list[fname]["field_type"] == "derived_func":
            # Read data
            try:
                if nonlocal_data is None:
                    data_convert = self.libyt.derived_func(grid.id, fname)
                else:
                    data_convert = nonlocal_data[grid.id][fname]
            except Exception as err:
                mylog.error("%s: %s", type(err).__name__, str(err))
                mylog.error(
                    "Cannot get derived field data in grid [%s] on MPI rank [%d]."
                    % (grid.id, grid.MPI_rank)
                )
                raise RuntimeError("libyt didn't get the data successfully.")
        else:
            # Since we only supports "cell-centered", "face-centered", "derived_func" tags for now
            # Raise an error if enter this block.
            raise ValueError(
                f"libyt does not have field_type [ {field_list[fname]['field_type']} ]"
            )

        # Swap axes or not, then return
        if field_list[fname]["contiguous_in_x"] is True:
            if data_convert.ndim == 3:
                return data_convert.swapaxes(0, 2)
            elif data_convert.ndim == 2:
                return np.expand_dims(data_convert.swapaxes(0, 1), axis=2)
            else:
                return np.expand_dims(data_convert, axis=(1, 2))
        else:
            return data_convert
