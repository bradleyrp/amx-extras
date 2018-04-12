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
multiply(nx=state.nx,ny=state.ny,nz=state.q('nz',1))
write_top('system.top')
if state.q('minimize',False): minimize('system')
else: copy_file('system.gro','system-minimized.gro')
copy_file('system-minimized.gro','system-residues.gro',)
write_structure_by_chain(structure='system-residues',gro='system')
gmx('make_ndx',ndx='system-groups',structure='system',inpipe="keep 0\nq\n",log='make-ndx-system-groups')
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

#! actually works to restart the step! rm -rf s03-large state.json && cp state_2.json state.json && make prep dextran_model_building_melt_multiply && cp expt_2.json expt.json && python -B script_2.py