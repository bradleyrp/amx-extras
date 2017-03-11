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
prepped_files = {'gro':state.here+get_last_gmx_call('pdb2gmx')['flags']['-o']}
state.protein_prepared = prepped_files
