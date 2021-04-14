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
        self._field_dtype = "float64"

###     ghost_zones != 0 is not supported yet
#       self.my_slice = (slice(ghost_zones,-ghost_zones),
#                        slice(ghost_zones,-ghost_zones),
#                        slice(ghost_zones,-ghost_zones))


    def _read_particle_coords(self, chunks, ptf):
        pass
#       chunks = list(chunks)   # generator --> list
#       p_idx  = self.ds.index._particle_indices

#       # shortcuts
#       par_posx = self._group_particle["ParPosX"]
#       par_posy = self._group_particle["ParPosY"]
#       par_posz = self._group_particle["ParPosZ"]

#       # currently libyt does not support multiple particle types
#       assert( len(ptf) == 1 )
#       ptype = list( ptf.keys() )[0]

#       for chunk in chunks:
#           for g1, g2 in particle_sequences(chunk.objs):
#               start = p_idx[g1.id    ]
#               end   = p_idx[g2.id + 1]
#               x     = np.asarray( par_posx[start:end], dtype=self._field_dtype )
#               y     = np.asarray( par_posy[start:end], dtype=self._field_dtype )
#               z     = np.asarray( par_posz[start:end], dtype=self._field_dtype )
#               yield ptype, (x, y, z)


    def _read_particle_fields(self, chunks, ptf, selector):
        pass
#       chunks = list(chunks)   # generator --> list
#       p_idx  = self.ds.index._particle_indices

#       # shortcuts
#       par_posx = self._group_particle["ParPosX"]
#       par_posy = self._group_particle["ParPosY"]
#       par_posz = self._group_particle["ParPosZ"]

#       # currently libyt does not support multiple particle types
#       assert( len(ptf) == 1 )
#       ptype   = list( ptf.keys() )[0]
#       pfields = ptf[ptype]

#       for chunk in chunks:
#           for g1, g2 in particle_sequences(chunk.objs):
#               start = p_idx[g1.id    ]
#               end   = p_idx[g2.id + 1]
#               x     = np.asarray( par_posx[start:end], dtype=self._field_dtype )
#               y     = np.asarray( par_posy[start:end], dtype=self._field_dtype )
#               z     = np.asarray( par_posz[start:end], dtype=self._field_dtype )

#               mask = selector.select_points(x, y, z, 0.0)
#               if mask is None: continue

#               for field in pfields:
#                   data = self._group_particle[field][start:end]
#                   yield (ptype, field), data[mask]


    def _read_chunk_data(self, chunk, fields):
        pass
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

        # TODO: Make this function return data-type that yt needs, modified if needed.
        mylog.debug("self.param_yt['field_list'] = %s", self.param_yt["field_list"])

        rv     = {}
        chunks = list(chunks)

        if selector.__class__.__name__ == "GridSelector":
            if not ( len(chunks) == len(chunks[0].objs) == 1 ):
                raise RuntimeError
            g = chunks[0].objs[0]
            for ftype, fname in fields:
                rv[(ftype, fname)] = self.grid_data[g.id][fname].swapaxes(0, 2)

            mylog.debug("###### (class IOHandlerlibyt, def _read_fluid_selection)")

            return rv

        if size is None:
            size = sum( (g.count(selector) for chunk in chunks for g in chunk.objs) )

        for field in fields:
            ftype, fname = field
            fsize = size
            rv[field] = np.empty(fsize, dtype=self._field_dtype)

        ng = sum( len(c.objs) for c in chunks )
        mylog.debug( "Reading %s cells of %s fields in %s grids",
                     size, [f2 for f1, f2 in fields], ng )

        offset = 0
        for chunk in chunks:
            for g in chunk.objs:
                for field in fields:
                    ftype, fname = field
### for ghost_zones != 0
#                   data_view = self.grid_data[g.id][fname][self.my_slice].swapaxes(0,2)
                    data_view = self.grid_data[g.id][fname][:,:,:].swapaxes(0,2)
                    offset   += g.select(selector, data_view, rv[field], offset)
        assert(offset == fsize)

        mylog.debug("###### (class IOHandlerlibyt, def _read_fluid_selection)")

        return rv

