"""
libyt-specific fields



"""

# -----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from yt.fields.field_info_container import FieldInfoContainer


# this is the FieldInfo subclass adopted when libyt fails to find a matched code frontend
class libytFieldInfo(FieldInfoContainer):
    # TODO[Polish]: Should append fields through libyt C-Library
    known_other_fields = ()
    known_particle_fields = ()

    def __init__(self, ds, field_list):
        super(libytFieldInfo, self).__init__(ds, field_list)

    def setup_fluid_fields(self):
        unit_system = self.ds.unit_system

    def setup_particle_fields(self, ptype):
        super(libytFieldInfo, self).setup_particle_fields(ptype)
