#!/usr/bin/env python

"""
MULTIPLY (a bilayer) to make it larger.
"""

from amx import *

init()
make_step(settings.step)
write_continue_script()
get_last_sources()
get_last_frame()
get_last_mdps()
multiply(nx=state.nx,ny=state.ny)
write_top('system.top')
bilayer_sorter(structure='system',ndx='system-groups')
if state.q('minimize',False): minimize('system')
else: copy_file('system.gro','system-minimized.gro')
equilibrate(structure='system-minimized',groups='system-groups')
