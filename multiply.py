#!/usr/bin/env python

"""
MULTIPLY

Function required for multiply scripts in automacs.
Currently designed for multiplying a bilayer system to make it larger.
"""

import re,os,shutil

def multiply(nx=1,ny=1,nz=1,quirky_ions=True,structure='system-previous'):
	"""
	Make a copy of a simulation box in multiple directions.
	Designed for membrane-(protein) systems.
	"""
	factor = nx*ny*nz
	#---collect composition details from the previous state settings object
	state.settings = state.before[-1]['settings']
	#---! previously loaded the prereqs directly into the state from the last settings
	state.lipids = state.before[-1]['lipids']
	state.force_field = state.before[-1]['force_field']
	shutil.rmtree(os.path.join(state.here,state.force_field)+'.ff')
	shutil.copytree(os.path.join(state.before[-1]['here'],state.force_field)+'.ff',
		os.path.join(state.here,state.force_field)+'.ff')
	settings['landscape_metadata'] = state.before[-1]['settings']['landscape_metadata']
	state.sol = state.before[-1]['settings']['sol']
	state.cation = state.before[-1]['settings']['cation']
	state.anion = state.before[-1]['settings']['anion']
	#---run genconf to multiply the system
	kwargs = {} if not state.buffer else {'flag':' -dist %.2f %.2f %.2f'%tuple(state.buffer)}
	#---we multiply the main system, but we may also need to multiply the restraints
	multiplies = [(structure,'system')]
	#---make sure the restraints were registered
	do_restraints = state.file_registry and 'system-leaflets-flat.gro' in state.file_registry
	#---add to the list of pairs of small-big files
	if do_restraints: multiplies += [('system-leaflets-flat','system-flat')]
	#---loop over configurations we need to multiply
	for from_fn,to_fn in multiplies:
		gmx('genconf',structure=from_fn,gro=from_fn+'-multiply-unordered',
			nbox="%d %d %d"%(nx,ny,nz),log='genconf-multiply-'+from_fn,**kwargs)
		struct = GMXStructure(state.here+from_fn+'-multiply-unordered'+'.gro')
		struct.regroup()
		struct.write(state.here+'%s.gro'%to_fn)
	#---modify the restraint rule
	if do_restraints:
		#---locate the restraint rule
		rule_ind = [ii for ii,i in enumerate(state.gmx_call_rules) 
			if i['command']=='grompp' and i['flag']=='r']
		if len(rule_ind)!=1: 
			raise Exception('there appear to be multiple gmx_call_rules for grompp restraints')
		else: rule = state.gmx_call_rules.pop(rule_ind[0])
		#---we update the restraint rule to use the modified bilayer
		rule['value'] = 'system-flat.gro'
		register_gmx_call(**rule)
		register_file('system-flat.gro')
	#---update the composition
	state.composition = struct.detect_composition()
	#---we retain the original lipid list since detect_composition cannot understand alternate names
	state.lipids = state.before[-1]['lipids']
