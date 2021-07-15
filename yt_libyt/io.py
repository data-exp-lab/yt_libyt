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
from itertools import groupby

from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.logger import ytLogger as mylog
from yt.geometry.selection_routines import AlwaysSelector



# group grids with consecutive indices together to improve the I/O performance
# --> grids are assumed to be sorted into ascending numerical order already
#def grid_sequences(grids):
#    for k, g in groupby( enumerate(grids), lambda i_x:i_x[0]-i_x[1].id ):
#        seq = list(v[1] for v in g)
#        yield seq

#def particle_sequences(grids):
#    for k, g in groupby( enumerate(grids), lambda i_x:i_x[0]-i_x[1].id ):
#        seq = list(v[1] for v in g)
#        yield seq[0], seq[-1]

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

        mylog.debug("self.grid_data.keys() = %s", self.grid_data.keys())

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

                    # g.id ptype particle number is 0 or g.id is not in local grids
                    # libyt.get_attr will return None
                    if x is None or y is None or z is None:
                        continue
                    else:
                        yield ptype, (x, y, z)


    def _read_particle_fields(self, chunks, ptf, selector):
        chunks = list(chunks)
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

                    # g.id ptype particle number is 0 or g.id is not in local grids
                    # libyt.get_attr will return None
                    if x is None or y is None or z is None:
                        continue

                    mask = selector.select_points(x, y, z, 0.0)
                    if mask is None:
                        continue

                    for field in ptf[ptype]:
                        data = self.libyt.get_attr(g.id, ptype, field)
                        # if ptype particle num in grid g.id = 0, get_attr will return None.
                        # then we shall continue the loop
                        # g.id ptype particle number is 0 or g.id is not in local grids
                        # libyt.get_attr will return None
                        if data is None:
                            continue
                        else:
                            yield (ptype, field), data[mask]


    def _read_chunk_data(self, chunk, fields):
        pass
        # TODO: Check this
#       rv = {}
#       if len(chunk.objs) == 0: return rv

#       for g in chunk.objs: rv[g.id] = {}

#       # Split into particles and non-particles
#       fluid_fields, particle_fields = [], []
#       for ftype, fname in fields:
#           if ftype in self.ds.particle_types:
#               particle_fields.append( (ftype, fname) )
#           else:
#               fluid_fields.append( (ftype, fname) )

#       # particles
#       if len(particle_fields) > 0:
#           selector = AlwaysSelector(self.ds)
#           rv.update( self._read_particle_selection(
#               [chunk], selector, particle_fields) )

#       # fluid
#       if len(fluid_fields) == 0: return rv

#       for field in fluid_fields:
#           ds = self._group_grid[ field[1] ]

#           for gs in grid_sequences(chunk.objs):
#               start = gs[ 0].id
#               end   = gs[-1].id + 1
#               data  = ds[start:end,:,:,:].transpose()
#               for i, g in enumerate(gs):
#                   rv[g.id][field] = np.asarray( data[...,i], dtype=self._field_dtype )
#       return rv


    def _read_fluid_selection(self, chunks, selector, fields, size):

        mylog.debug("#FLAG#")
        mylog.debug("yt/frontends/libyt/io.py (class IOHandlerlibyt, def _read_fluid_selection)")
        mylog.debug("fields = %s", fields)

        mylog.debug("self.param_yt['field_list'] = %s", self.param_yt["field_list"])
        field_list = self.param_yt["field_list"]
        rv     = {}
        chunks = list(chunks)

        # TODO: Need careful check for this if block
        if selector.__class__.__name__ == "GridSelector":
            if not ( len(chunks) == len(chunks[0].objs) == 1 ):
                raise RuntimeError
            g = chunks[0].objs[0]
            for ftype, fname in fields:
                rv[(ftype, fname)] = self.grid_data[g.id][fname].swapaxes(0, 2)

            mylog.debug("###### (class IOHandlerlibyt, def _read_fluid_selection)")

            return rv

        # TODO: Need careful check for this if block
        if size is None:
            size = sum( (g.count(selector) for chunk in chunks for g in chunk.objs) )

        for field in fields:
            ftype, fname = field
            fsize = size
            rv[field] = np.empty(fsize, dtype=self._field_dtype)

        ng = sum( len(c.objs) for c in chunks )
        mylog.debug( "Reading %s cells of %s fields in %s grids",
                     size, [f2 for f1, f2 in fields], ng )

        for field in fields:
            offset = 0
            ftype, fname = field
            mylog.debug("ftype, fname = %s", field)
            for chunk in chunks:
                for g in chunk.objs:
### for ghost_zones != 0
#                   data_view = self.grid_data[g.id][fname][self.my_slice].swapaxes(0,2)
                    # TODO: self.grid_data has all the g.id as keys, so we probably need additional check to prevent not parallel
                    if field_list[fname]["field_define_type"] == "cell-centered":
                        mylog.debug("self.grid_data[g.id][fname].shape = %s", self.grid_data[g.id][fname].shape)
                        data_convert = self.grid_data[g.id][fname][:, :, :]
                    elif field_list[fname]["field_define_type"] == "face-centered":
                        # TODO: This function is not quite correct, if size of patch is not square.
                        # convert to cell-centered
                        data_temp = self.grid_data[g.id][fname]
                        axis = np.argmax(data_temp.shape)
                        if axis == 0:
                            data_convert = 0.5 * (data_temp[:-1,:,:] + data_temp[1:,:,:])
                        if axis == 1:
                            data_convert = 0.5 * (data_temp[:,:-1,:] + data_temp[:,1:,:])
                        if axis == 2:
                            data_convert = 0.5 * (data_temp[:,:,:-1] + data_temp[:,:,1:])
                    elif field_list[fname]["field_define_type"] == "derived_func":
                        data_convert = self.libyt.derived_func(g.id, fname)
                    else:
                        # Since we only supports "cell-centered", "face-centered", "derived_func" tags for now
                        # Raise an error if enter this block.
                        raise ValueError("libyt does not have field_define_type [ %s ]" %
                                         (field_list[fname]["field_define_type"]))

                    # Swap axes or not
                    if field_list[fname]["swap_axes"] == True:
                        data_view = data_convert.swapaxes(0,2)

                    offset   += g.select(selector, data_view, rv[field], offset)
        assert(offset == fsize)

        mylog.debug("###### (class IOHandlerlibyt, def _read_fluid_selection)")

        return rv

