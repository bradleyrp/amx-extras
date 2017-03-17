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
prepped_files = {'gro':state.here+'vacuum-minimized.gro'}
state.protein_prepared = prepped_files
