#!/usr/bin/env python

"""
RESTART a GROMACS simulation from minimal components.
"""

from amx import *
import glob,re,shutil,os,stat

def unique_input(suffix,where='inputs'):
	"""Check for a unique file by suffix."""
	has_suffix = r'^(.*)\.%s$'
	fns = glob.glob(os.path.join(where,'*'))
	fns = [os.path.realpath(i) for i in fns 
		if re.match(has_suffix%suffix,os.path.basename(i))]
	if len(fns)==1: return fns[0]
	else: return None

make_step(settings.step)

# detect a single cpt and tpr in inputs
if settings.detect=='inputs':
	cpt = unique_input('cpt',where=settings.detect)
	tpr = unique_input('tpr',where=settings.detect)
	if not cpt or not tpr:
		raise Exception('cannot find a CPT and TPR file in %s'%settings.detect)
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

#! sidestep the error in which the cpt step is too far
if time_start!=0 and instruct_extend=='extend' and not settings.hard:
	instruct_extend = 'until'
	when = when + time_start

# standard restart requires a cpt/tpr only
tpr_start = 'md.part0001.tpr'
tpr_start = 'md.part0001.tpr'
mdp_start = 'md.part0001.mdp'
checkpoint_in = 'start.cpt'
mdrun_kwargs = {}
if not settings.hard:
	# convert the TPR
	convert_kwargs = {instruct_extend:when}
	gmx('convert-tpr',input=tpr,
		log='continue-tpr',out=tpr_start,**convert_kwargs)
	shutil.copyfile(cpt,state.here+checkpoint_in)
	mdrun_kwargs['checkpoint_in'] = checkpoint_in
# hard restart from a derived structure with new mdp options
else:
	# detects: mdp, top, cpt, tpr and optional ndx and ff folders
	#! dev: consider a json file to annotate the inputs
	mdp = unique_input('mdp',where=settings.detect)
	top = unique_input('top',where=settings.detect)
	if not mdp or not top: 
		raise exception('we require an mdp and top at %s'%settings.detect)
	ndx = unique_input('ndx',where=settings.detect)
	# collect force field folders
	ffs = [os.path.realpath(fn) 
		for fn in glob.glob(os.path.join(settings.detect,'*'))
		if re.match(r'^.+\.ff$',fn)
		and os.path.isdir(fn)]
	links = ffs+[mdp,top]
	if ndx: links += [ndx]
	for fn in links: 
		os.symlink(fn,os.path.join(state.here,os.path.basename(fn)))
	# built a structure from the checkpoint
	structure = 'start.gro'
	kwargs = dict(inpipe='0\n')
	gmx('trjconv',structure=cpt,input=tpr,
		out=structure,log='trjconv-cpt-structure',**kwargs)
	kwargs = dict(out=tpr_start,structure=structure,parameters=mdp,
		parameters_out=mdp_start,topology=top)
	if ndx: kwargs['groups'] = ndx
	gmx('grompp',log='grompp-restart',**kwargs)

# provide a continue script with a default extension time
fn_continue = state.here+'continue.sh'
with open('amx/gromacs/continue_std.sh') as fp: text = fp.read()
# substitute instructions so that the continuation script
#   either continues to an "until" mark or repeats the "extend" from the expt
text = re.sub('^MODE=(?:.+)$',text,
	'MODE=%s'%instruct_extend.upper(),flags=re.M)
if instruct_extend=='extend':
	text = re.sub('^EXTEND=(?:.+)$',text,
		'MODE=%s'%str(settings.extend),flags=re.M)
elif instruct_extend=='until':
	text = re.sub('^UNTIL=(?:.+)$',text,
		'MODE=%s'%str(settings.until),flags=re.M)
else: raise Exception('extension instructions must be "until" or "extend"')
with open(fn_continue,'w') as fp: fp.write(text)
st = os.stat(fn_continue)
os.chmod(fn_continue,st.st_mode | stat.S_IEXEC)

gmx('mdrun',
	input=tpr_start,
	log='mdrun-part0001',
	# deffnm ensures separate checkpoints
	deffnm='md.part0001',
	**mdrun_kwargs)
