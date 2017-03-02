#!/usr/bin/env python

"""
STRUCTURE EXPLORER
make cool new structures?
"""

from amx import *
init()
make_step(settings.step)
write_mdp()
make_coiled_coil(gro='structure-crude')
minimize('vacuum')
solvate('vacuum-minimized',gro='solvate',edges=state.water_edges,center=True)
write_top('solvate.top')
minimize('solvate')
counterions(structure='solvate',gro='counterions')
counterion_renamer('counterions')
minimize('counterions')
copy_file('counterions-minimized.gro','system.gro')
write_top('system.top')
grouper(protein=True,lipids=False)
equilibrate(groups='system-groups.ndx')
