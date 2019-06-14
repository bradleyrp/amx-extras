#!/usr/bin/env python

"""
RESTART a GROMACS simulation from minimal components.
"""

from amx import *

make_step(settings.step)

import glob,re,shutil
# detect a single cpt and tpr in inputs
if settings.detect=='inputs':
	has_suffix = r'^(.*)\.%s$'
	fns = glob.glob('inputs/*')
	fns_tpr,fns_cpt = [[os.path.realpath(i) for i in fns 
		if re.match(has_suffix%suf,os.path.basename(i))]
		for suf in ['tpr','cpt']]
	if len(fns_tpr)!=1 or len(fns_cpt)!=1:
		raise Exception('cannot find unique cpt/tpr: %s, %s'%(
			list(fns_cpt),list(fns_tpr)))
	cpt,tpr = fns_cpt[0],fns_tpr[0]	
else: raise Exception('unclear settings.detect: %s'%settings.detect)

ans = gmx('dump',cp=cpt)
time_start = re.findall(r'^t\s*=\s*(.*?)$',
	ans['stdout'],flags=re.M+re.DOTALL)
if len(time_start)!=1: 
	raise Exception('could not extract time from CPT dump')
time_start = float(time_start[0])

# handle continuation logic
instruct_extend = {}
if sum([settings.until!=None,settings.extend!=None])!=1: 
	raise Exception(
		'use only "extend" or "until" in the settings: %s'%settings)
if settings.until:
	instruct_extend = 'until'
	when = settings.until
elif settings.extend:
	instruct_extend = 'extend'
	when = settings.extend
else: raise Exception('dev')

"""
note that GROMACS throws an error if the cpt file step number is too late"
	The input requested 20000 steps, however the checkpoint file has already
	reached step 19498950. The simulation will not proceed, because either your
	simulation is already complete, or your combination of input files don't
	match.
to manage this, you should either start the simulation from a new structure
and MDP file (ideally, using a soon-to-be developed)
""


if time_start!=0 and instruct_extend=='extend':
	instruct_extend = 'until'
	when = when + time_start

# convert the TPR
convert_kwargs = {instruct_extend:when}
checkpoint_in = 'start.cpt'
gmx('convert-tpr',input=tpr,
	log='continue-tpr',out='start.tpr',**convert_kwargs)
shutil.copyfile(cpt,state.here+checkpoint_in)

# run mdrun
if not settings.legacy_part_names:
	"""
	the new default part names are:
		traj.part0002.trr
		traj_comp.part0002.xtc
		ener.part0002.edr
		md.part0002.log
	and gromacs will autodetect as long as you use the call below
	"""
	gmx('mdrun',
		input='start.tpr',
		checkpoint_in=checkpoint_in,
		noappend=True,
		log='mdrun-part0002')
# legacy naming scheme for a restart
else:
	base = 'md.part0002'
	gmx('mdrun',noappend=True,
		checkpoint_in=checkpoint_in,
		compressed='%s.xtc'%name,
		energy='%s.edr'%name,
		logfile='%s.log'%name,
		structure='%s.gro'%name,
		trajectory='%s.trr'%name,
		log='mdrun-part0002')
