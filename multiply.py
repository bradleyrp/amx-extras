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
	state.force_field = state.before[-1].get('force_field',state.settings.get('force_field'))
	shutil.rmtree(os.path.join(state.here,state.force_field)+'.ff')
	shutil.copytree(os.path.join(state.before[-1]['here'],state.force_field)+'.ff',
		os.path.join(state.here,state.force_field)+'.ff')
	try: settings['landscape_metadata'] = state.before[-1]['settings']['landscape_metadata']
	except: pass
	try: state.sol = state.before[-1]['settings']['sol']
	except: pass
	try: state.cation = state.before[-1]['settings']['cation']
	except: pass
	try: state.anion = state.before[-1]['settings']['anion']
	except: pass
	for itp in state.q('itp',[]):
		shutil.copyfile(state.before[-1].here+itp,state.here+itp)
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
		gmx_register_call(**rule)
		register_file('system-flat.gro')
	#---update the composition
	#! for the experiment script multiply-general.py we have to pass through a correction to the composition
	state.composition = struct.detect_composition(composition_adjust=settings.get('composition_adjust',None))
	#! hack for multiply_general expt
	renamer = settings.get('rename_detected_composition',None)
	if renamer: state.composition = [[renamer.get(i,i),j] for i,j in state.composition]
	#---we retain the original lipid list since detect_composition cannot understand alternate names
	try: state.lipids = state.before[-1]['lipids']
	except: pass

#! extremely bad that you cannot import this. it might be worth centralizing too
def write_structure_by_chain(structure='system',gro='system_chains'):
	"""
	Ensure that each chain has its own residue number. Useful for visualization.
	"""
	import numpy as np
	struct = GMXStructure(state.here+'%s.gro'%structure)
	polymer_molname = settings.molecule_name
	residue_name = settings.residue_name
	#! previously used DMR in the line below and DEX below that, in the component on the if not
	inds = np.where(struct.residue_names==residue_name)[0]
	#! assume that we only have a uniform melt of n_p, three beads per, minus two from the end points
	n_p = state.n_p
	bpm = settings.beads_per_monomer
	#! removed a -2 correction to the n_p for use in maltoheptaose which was necessary for removing terminals
	if not float(len(inds))/float(bpm)/(n_p)==component(polymer_molname): 
		raise Exception('failed to write structure by chains')
	resnums = (np.floor(np.arange(len(inds))/((n_p)*bpm))+1).astype(int)
	struct.residue_indices[inds] = resnums
	struct.write(state.here+'%s.gro'%gro)
