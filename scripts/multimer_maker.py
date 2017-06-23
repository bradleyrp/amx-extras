#!/usr/bin/env python

"""
MULTIMER MAKER
one-stop shop for complicated proteins
"""

from amx import *
init()
make_step(settings.step)
#write_mdp()
ProteinMultimer()
#minimize('vacuum')
#os.system('vmd %s'%(state.here+'multimer.gro'))
