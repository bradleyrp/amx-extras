#!/usr/bin/env python

import os,sys
import numpy as np

#---! you have to include not_reported functions from other places too?
_not_reported = ['crude_martini_residues','make_coiled_coil_backbone_crick','review_beads','vecnorm']

def make_coiled_coil_backbone_crick(n_residues,residue_shift=0,offset=0.0,overrides={}):
	"""
	The "simple" recipe corresponds to a Crick-style coiled coil.
	"""
	sqrt,cos,sin,tan,arctan,pi = np.sqrt,np.cos,np.sin,np.tan,np.arctan,np.pi

	#---describe a coiled coil with units in Angstroms
	specs = {
		'pitch':143,
		'minor_radius':1.495,
		'major_radius':4.7,
		'pitch_alpha_helix':5.4,
		'residues_per_turn':3.6,}
	specs.update(**overrides)
	pitch = specs['pitch']
	r0,r1 = specs['major_radius'],specs['minor_radius']
	winding = pitch/specs['pitch_alpha_helix']
	###
	### !!! YOU ABSOLUTELY MUST CHECK THE HANDEDNESS. SWITCHED TO "LEFT" BELOW BUT THAT IS BACKWARDS !!!
	###
	def helix(t,r=r0,p=pitch,offset=0.0,left=True): 
		sign = -1. if left else 1.
		return np.array([r*cos(sign*(t+offset)),r*sin(sign*(t+offset)),p/(2.*pi)*t])
	def frenet(t,r=r0,p=pitch,left=True):
		c = r/sqrt(r**2+p**2)
		s = p/sqrt(r**2+p**2)
		sign = 1 if not left else -1
		frame = np.array([
			[-c*sin(sign*t),-cos(sign*t),s*sin(sign*t)],
			[c*cos(sign*t),-sin(sign*t),-s*cos(sign*t)],
			[s,0,c],
			])
		return frame
	def coiled_coil(t,r=r1,a=winding):
		vec = np.array([0,cos(a*t),sin(a*t)]).T
		return r1*np.dot(frenet(t),vec)
	#---handle residue offsets in the spacings since that is the most direct way
	#---! note that doing a two helix bundle with maximum residue shift does not appear flush?
	spacings = 2*pi/winding/specs['residues_per_turn']*np.arange(residue_shift,n_residues+residue_shift)
	backbone = np.array([coiled_coil(p)+helix(p,offset=offset) for p in spacings]).T
	#---outgoing conversion from Angstroms to nanometers
	return backbone/10.0

def crude_martini_residues(backbone,seq,offset=0.0):
	"""
	"""
	lib_res_dn = state.residue_library
	#---! systematize the path to the amino acids ITP using the landscape or force field meta or something
	model = GMXTopology(state.here+os.path.join(state.force_field+'.ff','martini-v2.2-aminoacids.itp'))
	#---categorize the residues
	beader_type = {}
	for mol,spec in model.molecules.items():
		#---if one atom then it's easy and probably alanine
		if len(spec['atoms'])==1: beader_type[mol] = 'one'
		#---two atoms
		elif len(spec['atoms'])==2: beader_type[mol] = 'two'
		#---available in the library
		elif os.path.isfile(os.path.join(lib_res_dn,'%s.gro'%mol)): 
			beader_type[mol] = 'several'
		else: raise Exception('cannot find it in the residue library?')
	missings = [i for i in model.molecules if i not in beader_type]
	if any(missings): raise Exception('missing %r from the residue list'%missings)
	#---load structures for the "complicated" residues
	reslib = dict([(i,GMXStructure(os.path.join(lib_res_dn,'%s.gro'%i)))
		for i in list(zip(*filter(lambda x:x[1]=='several',beader_type.items())))[0]])
	#---each backbone point is the first in a list of residues
	pts = [[b] for b in backbone.T]
	#---best-guess positions of each residue
	for jj in range(1,len(pts)-1):
		resname = seq[jj]
		#---determine starting direction for non-terminal residues
		backbone_triple = np.array([pts[jj+i][0] for i in [-1,0,1]])
		#---! check if the points are in a straight line and if so make the point_out random
		#---project the middle point onto the line between flanking points, take its opposite
		a,p,b = backbone_triple
		point_out = vecnorm(-1*(a+np.dot(p-a,b-a)/np.dot(b-a,b-a)*(b-a)-p))
		#---build the residue
		bonds = model.molecules[resname].get('bonds')
		atoms = model.molecules[resname].get('atoms')
		constraints = model.molecules[resname].get('constraints')
		#---simple residues are constructed here
		if beader_type[resname] == 'one': continue
		if beader_type[resname] == 'two': 
			bond_length = float(bonds[0]['length']) if bonds else 0.35
			pts[jj] = np.concatenate((pts[jj],[bond_length*point_out+pts[jj][0]]))
		#---complicated residues are rotated onto the pointing-outward vector
		else: 
			#---backbone bead
			base = pts[jj][0]
			#---minimized points from the residue library
			xyz = np.array(reslib[resname].points)
			#---derived from the place-on-linker routine in the polymer melt code
			#---move the residue to the origin
			xyz -= xyz[0]
			#---backbone-to-first bead vector
			resvec = xyz[1]-xyz[0]
			#---compute the rotation matrix to align the residue with point_out
			a,b = vecnorm(resvec),vecnorm(point_out)
			rotation_angle = np.arccos(np.dot(a,b))
			rotation_axis = np.cross(a,b)
			rotation = rotation_matrix(rotation_axis,rotation_angle)
			rotated = np.dot(rotation,xyz.T).T
			shifted = rotated + base
			#---also replace the backbone residue
			pts[jj] = shifted
	#---we use flanking pts for the termini so we drop them here
	pts = pts[1:-1]
	return pts

def review_beads(*strands):
	"""
	Note that scaling factors assume nanometers.
	"""
	import mayavi
	from mayavi import mlab
	def pointplot(x,y,z,colormap='Spectral',r=1.0,color=(0,0,0)):
		mlab.points3d(x,y,z,colormap='Spectral',scale_factor=r,color=color)
	def lineplot(x,y,z,colormap='Spectral',r=0.1,color=(0,0,0)):
		mlab.plot3d(x,y,z,tube_radius=r,color=color)
	#---loop over "strands" (distinct protein structures)
	for beads in strands:
		pts_flat = np.array([i for j in beads for i in j])
		pointplot(*pts_flat.T,r=0.1)
		lineplot(*np.array([b[0] for b in beads]).T,r=0.05)
	mlab.show()

def make_coiled_coil_deprecated(gro):
	"""
	Replaces lib_helices.place_helices from previous incarnations of automacs.
	"""
	aa_codes1,aa_codes3 = GMXTopology._aa_codes1,GMXTopology._aa_codes3
	seq = map(lambda x:dict(zip(aa_codes1,aa_codes3))[x],state.sequence)
	#---sub-slice the sequence
	if state.sequence_slice:
		seqcut = [int(d)-1 for d in state.sequence_slice.split('-')]
		#---convert to pythonic slice
		seqcut = slice(seqcut[0],seqcut[1]+1)
	else: seqcut = slice(None,None)
	seq = list(seq)[seqcut]
	n_residues = len(seq)

	anti = state.oligomer_specs.get('anti',False)
	seq_a,seq_b = list(seq),list(seq)[::-1 if anti else 1]
	if state.oligomer_specs['n_coils']!=2: raise Exception('n_coils must be 2')
	#---crude_martini_residues must remove termini so we add false termini here
	backbone_a = make_coiled_coil_backbone_crick(n_residues+2)
	strand_a = crude_martini_residues(backbone_a,['ALA']+seq_a+['ALA'])
	backbone_b = make_coiled_coil_backbone_crick(n_residues+2,offset=np.pi,
		residue_shift=state.oligomer_specs['residue_shift'])
	strand_b = crude_martini_residues(backbone_b,['ALA']+seq_b+['ALA'])
	seqs = {'a':seq_a,'b':seq_b}

	#---check your work in three dimensions
	if state.review_3d: review_beads(strand_a,strand_b)

	#---! systematize the path to the amino acids ITP using the landscape or force field meta or something
	model = GMXTopology(state.here+os.path.join(state.force_field+'.ff','martini-v2.2-aminoacids.itp'))
	#---write the helices
	for snum,(strand,letter) in enumerate(zip([strand_a,strand_b],'ab')):
		#---compile the necessary metadat about the protein beads
		residue_indices,residue_names,atom_names,coords = [],[],[],[]
		residue_indices = np.concatenate([np.ones(len(j))*jj for jj,j in enumerate(strand)])+1
		residue_names = [seqs[letter][jj] for jj,j in enumerate(strand) for i in j]
		#---these are the atom names from the topology but MARTINI uses a simpler BB,SCN format
		atom_names = [i['atom'] for s in seqs[letter] for i in model.molecules[s]['atoms']]
		atom_names = [{0:'BB',1:'SC1',2:'SC2',3:'SC3',4:'SC4',5:'SC5'}[ii] 
			for s in seqs[letter] for ii,i in enumerate(model.molecules[s]['atoms'])]
		coords = np.array([i for j in strand for i in j])
		#---no negative coordinates
		coords -= coords.min(axis=0)
		#---box vectors slightly larger
		box_vector = coords.max(axis=0)*1.5
		#---avoid small-box errors
		box_vector[np.where(box_vector<10)] = 10.0
		#---use GMXStructure to write the gro
		another = GMXStructure(pts=coords+box_vector/2.0,box=box_vector,
			residue_indices=residue_indices,residue_names=residue_names,atom_names=atom_names)
		if snum==0: struct = another
		else: struct.add(another)
		another.write(state.here+'helix_%s.gro'%letter)
	struct.write(state.here+'helices.gro')
	#---! everything picks up in another function which is a fairly crude way to make the helix
	#---send the no-suffix structure names to the next step
	make_coiled_coil_simple('helix_a','helix_b')

def make_coiled_coil_SEE_MULTIMER(gro):
	"""
	Replaces lib_helices.place_helices from previous incarnations of automacs.
	Original version lacked Multimer -- this was developed to replace it with a more coherent protein class.
	"""
	aa_codes1,aa_codes3 = GMXTopology._aa_codes1,GMXTopology._aa_codes3
	seq = map(lambda x:dict(zip(aa_codes1,aa_codes3))[x],state.sequence)
	#---sub-slice the sequence
	if state.sequence_slice:
		seqcut = [int(d)-1 for d in state.sequence_slice.split('-')]
		#---convert to pythonic slice
		seqcut = slice(seqcut[0],seqcut[1]+1)
	else: seqcut = slice(None,None)
	seq = list(seq)[seqcut]
	n_residues = len(seq)

	anti = state.oligomer_specs.get('anti',False)
	seq_a,seq_b = list(seq),list(seq)[::-1 if anti else 1]
	if state.oligomer_specs['n_coils']!=2: raise Exception('n_coils must be 2')
	#---crude_martini_residues must remove termini so we add false termini here
	backbone_a = make_coiled_coil_backbone_crick(n_residues+2)
	strand_a = crude_martini_residues(backbone_a,['ALA']+seq_a+['ALA'])
	backbone_b = make_coiled_coil_backbone_crick(n_residues+2,offset=np.pi,
		residue_shift=state.oligomer_specs['residue_shift'])
	strand_b = crude_martini_residues(backbone_b,['ALA']+seq_b+['ALA'])
	seqs = {'a':seq_a,'b':seq_b}

	#---check your work in three dimensions
	if state.review_3d: review_beads(strand_a,strand_b)

	#---! systematize the path to the amino acids ITP using the landscape or force field meta or something
	model = GMXTopology(state.here+os.path.join(state.force_field+'.ff','martini-v2.2-aminoacids.itp'))
	#---write the helices
	for snum,(strand,letter) in enumerate(zip([strand_a,strand_b],'ab')):
		#---compile the necessary metadat about the protein beads
		residue_indices,residue_names,atom_names,coords = [],[],[],[]
		residue_indices = np.concatenate([np.ones(len(j))*jj for jj,j in enumerate(strand)])+1
		residue_names = [seqs[letter][jj] for jj,j in enumerate(strand) for i in j]
		#---these are the atom names from the topology but MARTINI uses a simpler BB,SCN format
		atom_names = [i['atom'] for s in seqs[letter] for i in model.molecules[s]['atoms']]
		atom_names = [{0:'BB',1:'SC1',2:'SC2',3:'SC3',4:'SC4',5:'SC5'}[ii] 
			for s in seqs[letter] for ii,i in enumerate(model.molecules[s]['atoms'])]
		coords = np.array([i for j in strand for i in j])
		#---no negative coordinates
		coords -= coords.min(axis=0)
		#---box vectors slightly larger
		box_vector = coords.max(axis=0)*1.5
		#---avoid small-box errors
		box_vector[np.where(box_vector<10)] = 10.0
		#---use GMXStructure to write the gro
		another = GMXStructure(pts=coords+box_vector/2.0,box=box_vector,
			residue_indices=residue_indices,residue_names=residue_names,atom_names=atom_names)
		if snum==0: struct = another
		else: struct.add(another)
		another.write(state.here+'helix_%s.gro'%letter)
	struct.write(state.here+'helices.gro')
	#---! everything picks up in another function which is a fairly crude way to make the helix
	#---send the no-suffix structure names to the next step
	make_coiled_coil_simple('helix_a','helix_b')

def backmapper_wrapper(target,destination,cwd,jitter=True,jitter_interval=None):
	"""
	Run the backmapper.
	!It would be awesome to rewrite the backmapper.
	!Backmapper code is python 2 only.
	"""
	target_fn = os.path.abspath(state.here+target)
	backward_path = os.path.join(os.getcwd(),state.backwards_script)
	#---! paths are absolute so we can execute backward.py in-place
	out = os.path.abspath(state.here+'backmap-'+destination)
	if sys.version_info>=(3,0): raise Exception('backward.py needs python 2')
	cmd = 'python %s -f %s -o %s.gro'%(backward_path,target_fn,out)
	print('[RUN] %s'%cmd)
	bash(cmd,cwd=cwd,log='backmapper-%s'%target)
	log_fn = os.path.join(state.here,'log-backmapper-%s'%target)
	with open(log_fn) as fp: logfile = fp.read()
	import re
	if re.search('(B|b)ailing',logfile):
		raise Exception('error in %s'%log_fn)
	if jitter: jitter_positions(structure=out,gro=destination+'.gro',
		cwd=state.here,interval=jitter_interval)

def jitter_positions(structure,gro,cwd='.',interval=None,nanjittter=0.1):

	"""
	The backmapper makes colinear points which chokes the minimizer until it segfaults.
	Here we add noise to the points to jitter them.
	"""

	if interval==None: interval = (-0.02,0.02)
	mol = GMXStructure(os.path.join(cwd,'%s.gro'%structure))
	mol.points += np.random.random(mol.points.shape)*float(interval[1]-interval[0])+interval[0]
	#---correct for any nan
	if not np.any(np.isnan(mol.points)): 
		mol.write(os.path.join(cwd,gro))
		return
	is_invalids = np.any(np.isnan(mol.points),axis=1)
	invalids = np.where(is_invalids)[0]
	#---take the nearst non-nan point
	nears = np.array([[np.argmin(np.abs(i-j)) for i in np.where(~is_invalids)][0] for j in invalids])
	#---jitter those points by an Angstrom
	mol.points[invalids] = mol.points[nears] + np.random.random((len(nears),3))*2.0-1.0
	mol.write(os.path.join(cwd,gro))

def make_coiled_coil_simple(*chains):
	"""
	Obviosuly this is just expedient. The idea and most code comes from script-helix-interactor.py.
	Takes a list of individual chains from whomever calls it.
	"""
	#---! note that we start a massive hack here: take crude CG points and backmap each chain separately
	#---! ...then combine these into one atomistic model and run that through martinize.py to get a topology
	#---! ...then return to the crude structure and minimize. THERE MUST BE A BETTER WAY TO DO THIS!
	#---loop over incoming chains
	for chain in chains: 
		#---run the backmapper
		backmapper_wrapper(chain+'.gro',chain+'-aa-back')
		#---we have to run pdb2gmx to make the residues whole again
		#---! backmapper misses atoms sometimes
		#---note the nice use of the "-p" flag below modifies the standard pdb2gmx so it includes a top file
		gmx('pdb2gmx',base='vacuum-aa-whole-'+chain,structure=chain+'-aa-back.gro',gro=chain+'-aa-whole-over',
			p='vacuum-aa-'+chain,log='pdb2gmx-'+chain,water=state.water,ff=state.atomistic_force_field)
		#---GRO lacks chains and anyway it's much easier to get topologies for each part separately
		#---! previously checked that the secondary structure matched the sequence before procedding
		gmx('editconf',structure=chain+'-aa-whole-over',gro=chain+'-aa-whole',c=True,d='2.0',
			log='editconf-center-%s'%chain)
		bash('python %s -f %s -x %s.pdb -o %s.top %s'%(
			os.path.abspath(state.martinize_script),chain+'-aa-whole.gro',
			chain+'-martinized',chain+'-martinized',
			'' if not state.secondary_structure else ' -ss %s'%state.secondary_structure),
			cwd=state.here,log='martinize')
		#---get the name of the itp from the topology
		top = state.here+chain+'-martinized.top'
		chain_top = GMXTopology(top)
		#---martinize always includes "martini.itp"
		itps = [i for i in chain_top.includes if i!='martini.itp']
		if len(itps)!=1: 
			import ipdb;ipdb.set_trace()
			raise Exception('non-unique itp referenced in %s'%top)
		itp_fn = itps[0]
		#---rename the molecule in the topology
		itp = GMXTopology(state.here+itp_fn)
		mols = itp.molecules.keys()
		if len(mols)!=1:
			raise Exception('non-unique list of molecules in this output from martinize: %s'%mols)
		#---rename the chains even if they are redundant
		new_name = 'protein_%s'%chain
		itp.molecule_rename(mols[0],new_name)
		#---add position restraints to the first BB atom
		#---! note that we simply reverse the sequence to make the antiparallel oligomer 
		#---! ...so we have to restrain the last position
		mol_spec = itp.molecules[new_name]
		posres_custom = {'funct': '1','fcy':'0','ai':'1','fcx':'0','fcz':'0'}
		posres_all = [dict(posres_custom,ai=i+1) for i in range(len(mol_spec['atoms']))]
		ind_end = {'helix_a':0,'helix_b':-1 if state.oligomer_specs.get('anti',False) else 0}[chain]
		first_bb = [ii for ii,i in enumerate(mol_spec['atoms']) if i['atom']=='BB'][ind_end]
		#---! hard-coded constraint force
		for i in 'xy': posres_all[first_bb]['fc%s'%i] = 1000
		itp.molecules[new_name]['position_restraints'] = posres_all
		itp.write(state.here+'protein_%s.itp'%chain,overwrite=True)
		if not state.itp: state.itp = []
		state.itp.append('protein_%s.itp'%chain)
		component('protein_%s'%chain,count=1)
		#---martinize writes pdbs so we change to gro
		#---awesome function: editconf has GRO.gro programmed as the output key but very easy to override
		gmx('editconf',f=chain+'-martinized.pdb',o=chain+'-martinized.gro',
			log='editconf-convert-%s'%chain)
	#---now that we have itp files and gro files we can reassemble them into one structure and topology
	for cnum,chain in enumerate(chains):
		another = GMXStructure(state.here+chain+'-martinized.gro')
		if cnum==0: struct = another
		else: struct.add(another)
	struct.write(state.here+'vacuum.gro')
	write_topology('vacuum.top')

	#---???
	if False:
		#---note the nice use of the "-p" flag below modifies the standard pdb2gmx so it includes a top file
		gmx('pdb2gmx',base='vacuum-aa',structure=chain+'.gro',gro='vacuum-aa-straight',
			p='vacuum-aa',log='pdb2gmx',water=state.water,ff=state.atomistic_force_field)
		if state.atomistic_minimization:
			raise Exception('dev')
			#---jitter whenever minimizing to avoid problems with colinear points
			lib_helices.jitter(structure='vacuum-aa-straight',gro='vacuum-aa',cwd=wordspace.step)
			minimize('vacuum-aa')
		else: copy_file('vacuum-aa-straight.gro','vacuum-aa-minimized.gro')
		#---! previously checked that the secondary structure matched the sequence before procedding
		bash('python %s -f %s.gro -x vacuum.pdb -o vacuum-martinize.top'%(
			os.path.abspath(state.martinize_script),'vacuum-aa-minimized')+
			('' if not state.secondary_structure else ' -ss %s'%state.secondary_structure),
			cwd=state.here,log='martinize')
	if False:
		wordspace.command_library['pdb2gmx']['-p'] = 'BASE.top'
		gmx('editconf',structure='helix',gro='helix-roomy',log='editconf-roomy',
			flag='-d %.1f'%wordspace.box_buffer)
		gmx('pdb2gmx',base='vacuum-aa',structure='helix-roomy.gro',gro='vacuum-aa-straight',top='vacuum-aa',
			log='pdb2gmx',water=wordspace['water'],ff=wordspace['atomistic_force_field'])
		if wordspace['atomistic_minimization']: 
			#---jitter whenever minimizing to avoid problems with colinear points
			lib_helices.jitter(structure='vacuum-aa-straight',gro='vacuum-aa',cwd=wordspace.step)
			minimize('vacuum-aa')
		else: filecopy(wordspace.step+'vacuum-aa-straight.gro',wordspace.step+'vacuum-aa-minimized.gro')
		seqcut = slice(*[int(d)+dd for dd,d in enumerate(wordspace.sequence_slice.split('-'))])
		assert len(wordspace.sequence[seqcut])==len(wordspace.secondary_structure),'ss must match sequence'
		bash('python %s -f %s.gro -x vacuum.pdb -o vacuum-martinize.top'%(
			os.path.expanduser(wordspace.martinize_path),'vacuum-aa-minimized')+
			('' if not wordspace['secondary_structure'] else ' -ss %s'%wordspace.secondary_structure),
			cwd=wordspace.step,log='martinize')
		gmx_run(gmxpaths['editconf']+' -f %s -o %s -d %.1f'%(
			'vacuum.pdb','vacuum.gro',wordspace.box_buffer),log='editconf-martinize')
	"""
	autodetect_protein_itp()
	component('Protein_A',count=2)
	write_top('vacuum.top')
	wordspace.mdp_specs['group'] = 'cgmd'
	write_mdp()
	minimize('vacuum')
	lib_helices.solvate('vacuum-minimized',gro='solvate')
	write_top('solvate.top')
	minimize('solvate')
	counterions(structure='solvate',gro='counterions')
	counterion_renamer('counterions')
	minimize('counterions')
	filecopy(wordspace.step+'counterions-minimized.gro',wordspace.step+'system.gro')
	write_top('system.top')
	equilibrate()
	if wordspace['run_part_two']=='yes': continuation() 
	else: write_continue_script()
	"""
	#import ipdb;ipdb.set_trace()

def place_helices_deprecated(gro,review=False,**kwargs):

	"""
	Draw a curve in 3-space and overlay an amino acid sequence on it.
	"""

	#---where to find minimized residue structures
	lib_res_dn = os.path.join(state.residue_library,'')
	#---uniform spacing of the backbone
	spacing = state.spacing

	#---! remove
	if 1:
		spacecurve = lambda x: (0,x**2/5.,x)
		spacecurve = lambda x: (x*np.cos(x),x*np.sin(x),0)
		spacecurve = lambda x: (1*np.cos(2*x),x,0)
		pitch,yaw = 0.5,4
		#---the following works
		spacecurve = lambda x,pitch=0.05,yaw=4:(pitch*np.cos(yaw*x),pitch*np.sin(yaw*x),x)
		spacecurve = lambda x,pitch=0.20,yaw=8.0:(pitch*np.cos(yaw*x),pitch*np.sin(yaw*x),x)
		#---broke: spacecurve = lambda x:np.array([x,0.0,0.0])
	else:
		#---define a spacecurve for the protein backbone
		#---get the spacecurve as a lambda function from the state
		if not type(state.spacecurve)==str: spacecurve = state.spacecurve
		else: spacecurve = eval(state.spacecurve)

	#---imports for 3D viewing
	if review:
		#---requires working mayavi installation (typically via qt4 backend)
		os.environ['ETS_TOOLKIT']='qt4'
		from mayavi import mlab
		import matplotlib as mpl
		import matplotlib.pylab as plt
		#---import locals
		from plotter3d import meshplot,meshpoints,pbcwire

	#---sub-slice the sequence
	if state.sequence_slice:
		seqcut = slice(*[int(d)-1 for d in state.sequence_slice.split('-')])
	else: seqcut = slice(None,None)

	#---interpret the "model" from an ITP file
	#model = CoarseGrainedModel('inputs/martini/martini.ff/martini-v2.2-aminoacids.itp')
	model = GMXTopology('inputs/martini/martini.ff/martini-v2.2-aminoacids.itp')
	
	#---get the 3-letter sequence
	seq = map(lambda x:dict(zip(model.aa_codes1,model.aa_codes3))[x],state.sequence)
	#---! TEMPORARY
	seq = seq[seqcut]
	#---make a list of residue points which we will infer from the model
	pts = [[[j*spacing for j in spacecurve(ii)]] for ii,i in enumerate(seq)]

	#---ROUTINE to draw the backbone as a 1D surface embedded in 3D
	#---independent parametric variable fed to spacecurve
	#---! need to automatically set this -- it should be larger than the sampled points
	#---! ...but this is non-trivial to predict so perhaps we let the user make sure
	t = np.arange(0,500,0.05)
	#---points along the spacecurve
	y = np.array([spacecurve(i) for i in t])
	#---distances between the points (note that these are probably not evenly-spaced)
	dist = np.sqrt(np.sum([(y[:-1,j]-y[1:,j])**2 for j in range(3)],axis=0))
	#---cumulative distance along the spacecurve
	dist_along = np.concatenate(([0],dist.cumsum()))
	#---use the scipy wrapper for FITPACK to create an interpolated spline
	spline,uuu = scipy.interpolate.splprep(y.T,u=dist_along,s=0)
	#---request points along the spline according to the sequence slice, at a particular spacing
	if not seqcut.start and not seqcut.stop: interp_d = np.arange(0,len(seq)+2)*spacing
	else: interp_d = np.arange(seqcut.start,seqcut.stop+2)*spacing
	#---sample the spline
	samps = np.transpose(scipy.interpolate.splev(interp_d,spline))
	if review:
		meshpoints(y,scale_factor=0.1,color=(1,1,1))
		meshpoints(samps,scale_factor=0.2)
	#---check that the spacings are close to the requested spacing
	spacings_actual = np.linalg.norm(samps[1:]-samps[:-1],axis=1)
	if spacings_actual.max()>spacing*5: 
		#---! better error reporting so the use knows to extend the input domain
		print('[ERROR] your input domain (t) is too small')
		trace()
	print('[NOTE] spacings: mean=%r,std=%r'%(spacings_actual.mean(),spacings_actual.std()))

	#---sometimes we want to reverse the direction (it's a long story)
	pts_both_directions = {}

	for backwards in [False,True]:

		#---use these sampled points as the backbone
		pts = [[s] for s in samps[::-1 if backwards else 1]]

		#---categorize the residues
		beader_type = {}
		for mol,spec in model.molecules.items():
			#---if one atom then it's easy and probably alanine
			if len(spec['atoms'])==1: beader_type[mol] = 'one'
			#---two atoms
			elif len(spec['atoms'])==2: beader_type[mol] = 'two'
			#---available in the library
			elif os.path.isfile(lib_res_dn+'%s.gro'%mol): beader_type[mol] = 'several'
		missings = [i for i in model.molecules if i not in beader_type]
		assert not missings,'missing %r from the residue list'%missing

		#---load structures for the "complicated" residues
		reslib = dict([(i,GMXStructure(lib_res_dn+'%s.gro'%i)) 
			for i in zip(*filter(lambda x:x[1]=='several',beader_type.items()))[0]])

		#---best-guess positions of each residue
		for jj in range(1,len(pts)-1):
			resname = seq[jj-1]
			#---determine starting direction for non-terminal residues
			backbone_triple = np.array([pts[jj+i][0] for i in [-1,0,1]])
			#---! check if the points are in a straight line and if so make the point_out random
			#---project the middle point onto the line between flanking points, take its opposite
			a,p,b = backbone_triple
			point_out = vecnorm(-1*(a+np.dot(p-a,b-a)/np.dot(b-a,b-a)*(b-a)-p))
			#---build the residue
			bonds = model.molecules[resname].get('bonds')
			atoms = model.molecules[resname].get('atoms')
			constraints = model.molecules[resname].get('constraints')
			#---simple residues are constructed here
			if beader_type[resname] == 'one': continue
			if beader_type[resname] == 'two': 
				bond_length = float(bonds[0]['length']) if bonds else 0.35
				pts[jj] = np.concatenate((pts[jj],[bond_length*point_out+pts[jj][0]]))
			#---complicated residues are rotated onto the pointing-outward vector
			else: 
				#---backbone bead
				base = pts[jj][0]
				#---minimized points from the residue library
				xyz = np.array(reslib[resname].points)
				#---derived from the place-on-linker routine in the polymer melt code
				#---move the residue to the origin
				xyz -= xyz[0]
				#---backbone-to-first bead vector
				resvec = xyz[1]-xyz[0]
				#---compute the rotation matrix to align the residue with point_out
				a,b = vecnorm(resvec),vecnorm(point_out)
				rotation_angle = np.arccos(np.dot(a,b))
				rotation_axis = np.cross(a,b)
				rotation = rotation_matrix(rotation_axis,rotation_angle)
				rotated = np.dot(rotation,xyz.T).T
				shifted = rotated + base
				#---also replace the backbone residue
				pts[jj] = shifted
		#---we use flanking pts for the termini so we drop them here
		pts = pts[1:-1]

		#---save these points
		pts_both_directions[{False:'fwds',True:'back'}[backwards]] = pts

	#---review the structure
	if review:
		#---plot
		meshpoints([i for j in pts for i in j],color=(1,1,1),scale_factor=0.2,opacity=0.5)

	structures_set = {}
	for direct in ['fwds','back']:
		pts = pts_both_directions[direct]
		#---compile the necessary metadat about the protein beads
		residue_indices,residue_names,atom_names,coords = [],[],[],[]
		residue_indices = np.concatenate([np.ones(len(j))*jj for jj,j in enumerate(pts)])+1
		residue_names = [seq[jj] for jj,j in enumerate(pts) for i in j]
		#---these are the atom names from the topology but MARTINI uses a simpler BB,SCN format
		atom_names = [i['atom'] for s in seq for i in model.molecules[s]['atoms']]
		atom_names = [{0:'BB',1:'SC1',2:'SC2',3:'SC3',4:'SC4',5:'SC5'}[ii] 
			for s in seq for ii,i in enumerate(model.molecules[s]['atoms'])]
		coords = np.array([i for j in pts for i in j])
		#---no negative coordinates
		coords -= coords.min(axis=0)
		#---box vectors slightly larger
		box_vector = coords.max(axis=0)*1.5
		#---avoid small-box errors
		box_vector[np.where(box_vector<10)] = 10.0

		lattice = kwargs.get('lattice',None)
		#---use GMXStructure to write the gro
		structures_set[direct] = GMXStructure(pts=coords+box_vector/2.0,box=box_vector,
			residue_indices=residue_indices,residue_names=residue_names,atom_names=atom_names)

	#---single placement
	if not lattice: structure.write(state.here+'%s.gro'%gro)

	from copy import deepcopy
	wind_backwards = kwargs.get('wind_backwards',False)
	direct = {True:'back',False:'fwds'}[wind_backwards[0]] if wind_backwards else 'fwds'
	structure = deepcopy(structures_set[direct])
	n_proteins = len(lattice)
	#---the first instance of the helix 
	structure.points += np.array(lattice[0])
	for nprot,xyz in enumerate(lattice[1:]):
		direct = {True:'back',False:'fwds'}[wind_backwards[1+nprot]] if wind_backwards else 'fwds'
		structure_add = deepcopy(structures_set[direct])
		for key in ['residue_indices','residue_names','atom_names']:
			structure.__dict__[key] = np.concatenate((structure.__dict__[key],structure_add.__dict__[key]))
		structure.points = np.concatenate((structure.points,structure_add.points + np.array(xyz)))
	structure.write(state.here+'%s.gro'%gro)
