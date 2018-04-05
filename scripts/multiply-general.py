#!/usr/bin/env python

"""
MULTIPLY a simulation to make it larger.
Forked from the bilayer multiply step.
"""

from amx import *

init()
make_step(settings.step)
write_continue_script()
get_last_sources()
gmx_get_last_frame()
get_last_mdps()
# adding mdp_specs to settings will overwrite them
if state.mdp_specs: write_mdp(param_file=state.mdp_parameters)
multiply(nx=state.nx,ny=state.ny)
write_top('system.top')
if state.q('minimize',False): minimize('system')
else: copy_file('system.gro','system-minimized.gro')
# the multiply procedure can benefit from a sequence of topology changes
if state.equilibrate_tops:
	# using a sequence of mdp/ff changes overrides the settings.equilibration
	for rnum,(mdp,top) in enumerate(state.equilibrate_tops):
		state.force_field = top
		write_topology('system.top')
		restart_clean(part=rnum+1,mdp='input-md-%s-eq-in'%mdp if mdp else 'input-md-in',
			structure='system-minimized' if rnum==0 else 'md.part%04d'%rnum,
			groups='system-groups')
# equilibration requires a sequence of MDP files from the previous step
else: equilibrate(structure='system-minimized',groups='system-groups')

#! iterative development restart with::: rm -rf s02-large state.json && cp state_1.json state.json && make prep exo70_control_pc_ps && cp expt_2.json expt.json && python -B script_2.py

#! this one: rm -rf s03-large state.json && cp state_1.json state.json && make prep dextran_model_building_melt && cp expt_2.json expt.json && python -B script_2.py
#! failed at this ... oops