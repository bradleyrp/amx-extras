#!/usr/bin/env python

"""
Multimer class.

Organize protein multimers.
Started development on 2017.5.3.
Development started for the purposes of organizing proteins for bilayer adhesion. 
May be expanded to other use-cases later?
"""

import os,sys,re,shutil,copy,glob
from ortho.dictionary import DotDict
import numpy as np
#---! this is required on some computers but not others possibly due to weird imports
from .geometry_tools import plane_project,vecangle

def reviewer3d(**kwargs):
	"""Debug 3D manipulations."""
	import mayavi
	from mayavi import mlab
	def pointplot(x,y,z,colormap='Spectral',r=1.0,color=(0,0,0)):
		mlab.points3d(x,y,z,colormap='Spectral',scale_factor=r,color=color)
	for key,val in kwargs.items():
		pointplot(*np.array(val['points']).T,r=val.get('r',1.0),color=val.get('color',(0,0,0)))
	mlab.show()

class PDBStructure:

	def __init__(self,source):
		"""
		Parse a PDB into a GMXStructure.
		"""
		self.source = source
		#---all PDBs are represented internally as GRO files
		self.basename = re.match('^(.*?)\.pdb$',source).group(1)
		#---! nice example of overrides in GMX below with the "-f" flag
		gmx('editconf',f=source,gro=self.basename,
			log='editconf-convert-%s'%self.basename)
		#---read coordinates with the AUTOMACS structure class
		self.gro = GMXStructure('%s%s.gro'%(state.here,self.basename))
		#---! parse the file for important details?
		#---! move to a particular location?

class ProteinMultimer:

	def __init__(self):
		"""
		Orchestrate the construction of a protein multimer for reception by other AUTOMACS components.
		"""
		#---primary input method is a list of dictionary commands from settings
		self.pseudoscript = settings.getcheck('pseudoscript',list)
		#---child objects are corralled in the contents
		self.contents,self.loci = {},{}
		self.main()
		#self.writehack()

	def writehack(self):
		"""TEMPORARY HACK."""
		#---! renumber is False to skip that step
		self.contents['body'].gro.write(state.here+'multimer.gro',renumber=False)

	def write_contents(self,**kwargs):
		"""
		Write a constituent structure.
		"""
		name = kwargs.pop('name')
		gro = kwargs.pop('gro')
		if kwargs: raise Exception('unprocessed kwargs: %s'%kwargs)
		subject = self.select_gro(name)
		subject.write(state.here+'%s.gro'%gro,renumber=False)

	def main(self,debug_pseudoscript=False):
		"""The main loop iterates over commands in the pseudoscript."""
		for cnum,(command,kwargs) in enumerate(self.pseudoscript): 
			if not debug_pseudoscript: getattr(self,command)(**kwargs)
			else:
				try: getattr(self,command)(**kwargs)
				except Exception as e:
					from runner.makeface import tracebacker
					raise Exception('failed in the pseudoscript, command number %d: %s,%s: %s'%(
						cnum+1,command,kwargs,e))

	def pdb(self,**kwargs):
		"""Load a PDB into the multimer."""
		name = kwargs.pop('name',None)
		source = kwargs.pop('source',None)
		if kwargs: raise Exception('unprocessed kwargs: %s'%kwargs)
		from runner.acme import get_path_to_module
		shutil.copyfile(get_path_to_module(source),state.here+'%s.pdb'%name)
		self.contents[name] = PDBStructure('%s.pdb'%name)

	def select(self,name,kind=None):
		"""Get an item from the contents with type checking."""
		if name not in self.contents: raise Exception('cannot find %s in contents'%name)
		classname = self.contents[name].__class__.__name__
		if kind and classname!=kind:
			Exception('found %s in contents but it has class %s and not %s'%(name,classname,kind))
		else: return self.contents[name]

	def select_gro(self,name):
		"""
		Return the wrapped GMXStructure from an object in contents.
		"""
		subject = self.select(name)
		if subject.__class__.__name__=='GMXStructure': struct = subject
		elif subject.__class__.__name__=='PDBStructure': struct = subject.gro
		else: raise Exception('need a GMXStructure from %s'%name)
		return struct

	def locus(self,**kwargs):
		"""
		Identify and save key position.
		!!! This has been superceded by a more expedient solution made with translate, etc.
		"""
		name = kwargs.pop('name',None)
		resid = int(kwargs.pop('resid',0))
		atom_name = kwargs.pop('atom_name',None)
		locus_name = kwargs.pop('locus_name',None)
		if not resid: raise Exception('locus needs a resid')
		if not atom_name: raise Exception('locus needs an atom_name')
		if not locus_name: raise Exception('locus needs a locus_name')
		if kwargs: raise Exception('unprocessed kwargs: %s'%kwargs)
		#---defining the locus depends on the type of the subject
		subject = self.contents.get(name)
		if not subject: raise Exception('cannot find %s in contents'%name)
		subject = self.select(name,kind='PDBStructure')
		matches = np.where(np.all((subject.gro.residue_indices==resid,
			subject.gro.atom_names==atom_name),axis=0))[0]
		if len(matches)==0: 
			raise Exception('failed to find a locus with resid %s and atom_name %s'%(resid,atom_name))
		elif len(matches)>1: 
			raise Exception('multiple points found for resid %s and atom_name %s'%(resid,atom_name))
		else: self.loci[locus_name] = {'name':name,'index':matches[0]}

	def lay_flat(self,**kwargs):
		"""
		Lay elongated coordinates "flat" along some direction.
		Requires @bilayers/codes/adhere_protein_bilayer.py in extensions 
		and lay_coords_flat in _shared_extensions there.
		"""
		name = kwargs.pop('name',None)
		direction = kwargs.pop('direction','y')
		if kwargs: raise Exception('unprocessed kwargs: %s'%kwargs)
		struct = self.select_gro(name)
		struct.points = lay_coords_flat(struct.points,direction=direction)

	def inspect_protein_rotations(self,protein,*args,**kwargs):
		"""
		Custom inspection for use at various points in the rotation.
		"""
		do_review = kwargs.pop('do_review',False)
		coords = protein.points
		out = dict(body=dict(points=coords,r=0.05))
		out.update(com=dict(points=[coords.mean(axis=0)],r=0.3,color=(1,1,1)))
		for aa,arg in enumerate(args):
			for dd,d in enumerate(arg): 
				out.update(**{str(dd):dict(points=protein.select_center(d),r=0.2,color=(1,0,0))})
			out.update(**{str('c%d'%aa):dict(points=protein.select_center(arg),r=0.2,color=(1,0,0))})
		for key,val in kwargs.items():
			out.update(key=dict(points=[val],r=0.2,color=(0,1,0)))
		if do_review: 
			reviewer3d(**out)
			sys.exit()

	def banana(self,**kwargs):
		"""
		Mimic almost verbatim the banana routine (minus lay_flat above).
		"""
		name = kwargs.pop('name',None)
		group_down = kwargs.pop('group_down')
		group_origin = kwargs.pop('group_origin',copy.deepcopy(group_down))
		direction = kwargs.pop('direction','y')
		direction_down = kwargs.pop('direction_down','z')
		if kwargs: raise Exception('unprocessed kwargs: %s'%kwargs)

		ref_axis,down_axis = np.zeros(3),np.zeros(3)
		ref_axis['xyz'.index(direction)] = 1
		down_axis['xyz'.index(direction_down)] = -1
		protein = self.select(name,kind='PDBStructure').gro
		#---align the principal axis with a particular direction
		self.lay_flat(name=name,direction=direction)
		#---begin transplant from place_protein_banana
		#---! try to remove redundancy when possible
		downer = protein.select_center(group_down)
		centroid = protein.cog(protein.select('all'))
		coords = np.array(protein.points)
		#---all rotations are about the mean
		coords -= centroid
		#---project the vector between centroid and downward-facing group onto the plane normal to the direction
		#---now that we are mean-centering, we do not want: axis = vecnorm(downer-centroid)
		axis = vecnorm(downer)
		#---! previously used: principal = vecnorm(principal_axis(coords))
		#---! ...however the lay flat procedure might not be perfect, so we use the reference axis instead
		principal = vecnorm(np.array(ref_axis))
		projected = vecnorm(plane_project(axis,principal))
		### self.inspect_protein_rotations(protein,group_down,projected=projected)
		#---rotate protein along the direction so the axis points down
		rotation_angle = vecangle(down_axis,projected)
		#---! removed a buggy multiplicative inverse on the rotation angle below
		rotation = rotation_matrix(principal,rotation_angle)
		#---apply the rotation
		coords_rotated = np.dot(coords,rotation)
		#---shift back to the original centroid now that rotation is complete
		coords_rotated += centroid
		#---write the modified structure
		protein.points = coords_rotated
		#self.inspect_protein_rotations(protein,group_down,projected=projected)
		#---require an origin group to act as the reference point for the protein
		#group_origin = state.q('group_origin')
		#if not group_origin: raise Exception('banana requires group_origin or group_down in the settings')
		#pts_origin = protein.select_center(group_origin)
		#protein.points -= pts_origin
		#---make sure the box is big enough if we want to check it in e.g. VMD
		protein.box = protein.points.ptp(axis=0)-protein.points.min(axis=0)+1.0
		#---! still need to apply "center_selection" (lipid pocket) and "over"
		#---! also need to recenter the protein laterally (see original code)
		#---! wait to write the GRO file
		#---! okay that the final data are mean-centered? apply the pocket lipid here?
		#---! debugging the rotation
		if False:
			reviewer3d(
				body=dict(points=coords,r=0.05),
				com=dict(points=[centroid],r=0.3,color=(1,1,1)),
				downer=dict(points=[downer],r=0.3,color=(1,0,1)),
				projected=dict(points=[projected+centroid],r=0.3,color=(1,1,0)),)
			downer_rotated = protein.select_center(group_down)
			axis_rotated = vecnorm(downer_rotated-centroid)
			principal_rotated = vecnorm(principal_axis(coords_rotated))
			projected_rotated = vecnorm(plane_project(axis_rotated,principal_rotated))
			reviewer3d(
				downer_rotated=dict(points=[downer_rotated],r=1.0,color=(1,0,1)),
				projected_rotated=dict(points=[projected_rotated+centroid],r=1.0,color=(1,1,0)),
				com=dict(points=[coords_rotated.mean(axis=0)],r=1.0,color=(1,1,1)),
				rotated=dict(points=coords_rotated,r=0.05,color=(1,1,1)),)
		if False:
			coords = protein.points
			kwargs = dict(body=dict(points=coords,r=0.05))
			kwargs.update(com=dict(points=[coords.mean(axis=0)],r=0.3,color=(1,1,1)))
			for dd,d in enumerate(group_down): 
				kwargs.update(**{str(dd):dict(points=protein.select_center(d),r=0.2,color=(1,0,0))})
			kwargs.update(down=dict(points=protein.select_center(group_down),r=0.2,color=(1,0,0)))
			reviewer3d(**kwargs)
			sys.exit()

	def make_coiled_coil(self,**kwargs):
		"""
		Make a coiled-coil object.
		"""
		name = kwargs.pop('name')
		sequence = kwargs.pop('sequence')
		oligomer_specs = kwargs.pop('oligomer_specs')
		if kwargs: raise Exception('unprocessed kwargs: %s'%kwargs)

		#---porting make_coiled_coil from structural_biology.py
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

		#---anti or parallel sequences
		anti = oligomer_specs.get('anti',False)
		seq_a,seq_b = list(seq),list(seq)[::-1 if anti else 1]
		if oligomer_specs['n_coils']!=2: raise Exception('n_coils must be 2')
		#---!!!!!!!!!!!!!!!
		overrides = {'major_radius':10.0}
		#---crude_martini_residues must remove termini so we add false termini here
		backbone_a = make_coiled_coil_backbone_crick(n_residues+2,overrides=overrides)
		strand_a = crude_martini_residues(backbone_a,['ALA']+seq_a+['ALA'])
		backbone_b = make_coiled_coil_backbone_crick(n_residues+2,offset=np.pi,
			residue_shift=oligomer_specs.get('residue_shift',0),overrides=overrides)
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
		#make_coiled_coil_simple('helix_a','helix_b')
		#---porting make_coiled_coil_simple from structural_biology.py
		chains = 'helix_a','helix_b'

		#---! note that we start a massive hack here: take crude CG points and backmap each chain separately
		#---! ...then combine these into one atomistic model and run that through martinize.py 
		#---! ...to get a topology
		#---! ...then return to the crude structure and minimize. THERE MUST BE A BETTER WAY TO DO THIS!
		#---loop over incoming chains
		for chain in chains: 
			#---run the backmapper
			backmapper_wrapper(chain+'.gro',chain+'-aa-back',cwd=state.here)
			#---we have to run pdb2gmx to make the residues whole again
			#---! backmapper misses atoms sometimes
			#---note the nice use of the "-p" flag below modifies the standard pdb2gmx so it includes 
			#---...a top file
			water_model = 'none' if state.water==None else state.water
			gmx('pdb2gmx',base='vacuum-aa-whole-'+chain,structure=chain+'-aa-back.gro',
				gro=chain+'-aa-whole-over',p='vacuum-aa-'+chain,log='pdb2gmx-'+chain,
				water=water_model,ff=state.atomistic_force_field)
			#---GRO lacks chains and anyway it's much easier to get topologies for each part separately
			#---! previously checked that the secondary structure matched the sequence before procedding
			gmx('editconf',structure=chain+'-aa-whole-over',gro=chain+'-aa-whole',c=True,d='2.0',
				log='editconf-center-%s'%chain)
			bash('python %s -f %s -x %s.pdb -o %s.top %s'%(
				os.path.abspath(state.martinize_script),chain+'-aa-whole.gro',
				chain+'-martinized',chain+'-martinized',
				'' if not state.secondary_structure else ' -ss %s'%state.secondary_structure),
				cwd=state.here,log='martinize-%s'%chain)
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
			ind_end = {'helix_a':0,'helix_b':-1 if oligomer_specs.get('anti',False) else 0}[chain]
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
		#---! renamed from vacuum to linker (the name) 
		#---! ... but this was probably used for vacuum simulation of the linker in a previous demo
		struct.write(state.here+'%s.gro'%name)
		write_topology('%s.top'%name)

		#---???
		if False:
			#---note the nice use of the "-p" flag below modifies the standard pdb2gmx 
			#---...so it includes a top file
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
			gmx('pdb2gmx',base='vacuum-aa',structure='helix-roomy.gro',
				gro='vacuum-aa-straight',top='vacuum-aa',log='pdb2gmx',
				water=wordspace['water'],ff=wordspace['atomistic_force_field'])
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

		#---when this is complete we register the result in the contents list
		#---! fix the vacuum.gro to use a better name
		self.contents[name] = GMXStructure(state.here+'%s.gro'%name)

	def callout(self,**kwargs):
		"""
		Interface the Multimer class with automacs functions. 
		This function packages a fake `state` variable for automacs. 
		This provides more granular control than the usual state.

        #('callout',{'function':'martinize','state':{
        #    'martinize_path':'inputs/martini/bin/martinize.py',
        #    'dssp_path':'inputs/martini/bin/dssp-2.0.4-linux-amd64',
        #    }}),
        #---! dropping the callout idea because it's too hard to recap the state
        #---! ...basic idea is to have the function accept the state as an argument otherwise
        #---! ...get the global state. this way it can be called as part of the automacs workflow
        #---! ...or as part of 

		"""
		state = DotDict(**kwargs.pop('state',{}))
		function = kwargs.pop('function')
		if kwargs: raise Exception('unprocessed kwargs: %s'%kwargs)
		#---! check for the function and check that it takes the state
		globals()[function](state=state)

	def martinize(self,**kwargs):
		"""
		Use martinize to generate a coarse-grained protein.
		Adapted and refined from martini/martini.py for calls by Multimer.
		"""
		structure = kwargs.pop('structure')
		gro = kwargs.pop('name')
		log = kwargs.pop('log','martinize')
		secondary_structure = kwargs.pop('secondary_structure',None)
		if kwargs: raise Exception('unprocessed kwargs %s'%kwargs)
		#---convert GRO to PDB because martinize likes it that way
		gmx('editconf',structure=structure,o='%s.pdb'%structure,log='editconf-martinize-%s'%structure)
		martinize_fn = state.martinize_path
		#---this function is run from the step but martinize_path is relative to root
		if not os.path.isfile(martinize_fn): raise Exception('cannot find martinize at %s'%martinize_fn)
		cmd = 'python '+os.path.abspath(martinize_fn)+' -v -p backbone '
		#---! write a separate position restrained version if desired by dropping -p None
		cmd += ' -f %s -o %s.top -x %s.pdb '%('%s.pdb'%structure,gro,gro)
		#---dssp use should be standard for MARTINI proteins
		if state.dssp_path and secondary_structure==None: 
			dssp_fn = os.path.abspath(os.path.expanduser(state.dssp_path))
			if not os.path.isfile(dssp_fn):
				raise Exception('cannot find %s'%dssp_fn)
			cmd += ' -dssp %s'%dssp_fn
		#---customized secondary structure
		elif secondary_structure!=None: cmd += ' -ss %s'%secondary_structure
		else: raise Exception('secondary structure problem!')
		if state.martinize_ff_version: cmd += ' -ff %s'%state.martinize_ff_version
		if state.martinize_flags: cmd += ' '+state.martinize_flags
		bash(cmd,cwd=state.here,log=log)
		if not os.path.isfile(state.here+'%s.pdb'%gro): raise Exception('martinize failed')
		gmx_run(state.gmxpaths['editconf']+' -f %s.pdb -o %s.gro'%(gro,gro),log='editconf-convert-pdb')
		#---end martinize sequence, truncated from the original, which loaded ITP files into state
		#---the martinize.py script says that "this should probably be gathered in a Universe class"
		#---...and we could not agree more
		#---parse the resulting topology file for ITP files
		#---note that martinize will not produce multiple redundant ITP files which is good
		#---...but we need to be careful not to do that either
		#---! protection against multiple redundant ITP files goes where?
		#---parse the resulting topology file for ITP files
		with open(state.here+'%s.top'%gro) as fp: topfile = fp.read()
		#---collect ITPs
		itps = re.findall('\#include\s+"(.*?)"',topfile,re.M+re.DOTALL)
		#---assume all ITP files on disk are due to martinize
		new_itps = [i for i in itps if os.path.isfile(state.here+i)]
		#---use the topology class to catalogue the files and their molecules
		assoc_topology = []
		for new_itp in new_itps:
			#---all ITPs built with martinize are renamed (prefixed with the output name and lowercased)
			topname = re.match('^(.+)\.itp',new_itp).group(1).lower()
			#---!!! VERY FRUSTRATED WITH INELEGANT DATA STRUCTURES SO HACKING! also martinize broken!
			with open(state.here+new_itp,'r') as fp: text = fp.read()
			deletions = [
				'#ifndef RUBBER_FC(.*?)#endif',
				'#ifndef NO_RUBBER_BANDS(.*?)#endif',
				'#ifndef POSRES_FC(.*?)#endif',
				'#ifdef POSRES(.*?)#endif',]
			for kill in deletions:
				text = re.sub(kill,'',text,flags=re.M+re.DOTALL)
			with open(state.here+new_itp,'w') as fp: fp.write(text)
			#---! different types implicitly refer to different files. this should be made explicit
			self.contents['%s_%s'%(gro,topname)] = GMXTopology(state.here+new_itp,
				defs={'posres':False,'no_rubber_bands':True},constraints_to_bonds=True)
			#---save a bonded version
			self.contents['%s_%s_bonded'%(gro,topname)] = GMXTopology(state.here+new_itp,
				defs={'posres':False,'no_rubber_bands':True},constraints_to_bonds=True)
			assoc_topology.append('%s_%s'%(gro,topname))
			#---to avoid clutter, and because martinize may be called repeatedly, we clean the files which
			#---...can easily be written later
			#---! later when the files are written we need to avoid name collisions in the moleculename column
			os.remove(state.here+new_itp)
		#---read the resulting structure back into the contents
		self.contents[gro] = GMXStructure(state.here+'%s.gro'%gro)
		#---let the user know which ITPs we found
		return assoc_topology

	def duplicate(self,name,new):
		"""
		"""
		self.contents[new] = copy.deepcopy(self.contents[name])

	def translate(self,**kwargs):
		"""
		!!! Under development from a monthslong gap.
		Since multimers are connected in multiple places we proceed with piecewise transformations
		"""
		if set(kwargs.keys())<=set([
			'reference_name','reference_selection','subject_name','subject_selection','offset']):
			#---move a reference to a subject
			ref_name = kwargs.pop('reference_name')
			ref_sel = kwargs.pop('reference_selection')
			subject_name = kwargs.pop('subject_name')
			subject_sel = kwargs.pop('subject_selection')
			ref = self.select(ref_name)
			subject = self.select(subject_name)
			#---the reference center is the position we are moving the subject selection to
			ref_center = ref.cog(ref.select(ref_sel))
			subject_center = subject.cog(subject.select(subject_sel))
			subject.points += ref_center - subject_center + kwargs.pop('offeset',[0.,0.,0.])
		else: raise Exception('invalid arguments')

	def rotate(self,**kwargs):
		"""
		Rotate some points
		"""
		ref_axis = np.array(kwargs.pop('axis',[0,0,1]))
		rotation_angle = kwargs.pop('angle')/180*np.pi
		subject = self.select(kwargs['name'])
		#---save the shift to return the object after rotating
		shift = subject.cog(subject.select('all'))
		subject.points -= shift
		rotation = rotation_matrix(ref_axis,rotation_angle)
		subject.points = np.dot(subject.points,rotation)
		#---return the points
		subject.points += shift

	def exo70_builder(self,**kwargs):
		"""
		PROTOTYPING CODE TO BUILD EXO70 MODELS!
		NEEDS TO BE GENERALIZED EVENTUALLY!

		Previously, we assembled two protein body elements and a linker.
		Everything is in the right position, so here we stitch it all together.
		Take the linker and split it. Write a structure with the body and then backmap.
		Use the backmapped structure to generate a topology.
		Martinize will read (and expunge the topology).
		We will internall track the link between structure and topology for now.
		"""
		n_mers = 2
		itp_names,itp_names_equilibrate = [],[]
		#---! for posterity
		state.protein_prepared = {'gro':None,'top':None,'composition':[]}
		#---feeling ambitious so this is going in a loop over mers
		for mono_num in range(n_mers):
			base = self.select('body_cg_%d'%(mono_num+1))
			link = self.select('linker')
			#---the linker includes chains for each monomer so we pull them out here
			#---! obviously this is somewhat arbitrary
			#---! change this order, the resids in the translate, and the rotation in the script
			#---! ... to change all of the 
			#---number of points in a linker monomer
			nplm = 191
			#---! INITIAL ORDER WAS WRONG AND THIS CAUSED UNTOLD PROBLEMS
			link_slice = [slice(201,None),slice(None,201)][mono_num]
			#---! ALSO HAD TO REVERSE IT BUT ONLY FOR THE FIRST ONE
			link_slice = [np.arange(nplm*1,nplm*2)[::-1],np.arange(nplm*0,nplm*1)[::1]][mono_num]
			#---concatenate the structures
			points = np.concatenate((link.points[link_slice],base.points))
			residue_names = np.concatenate((link.residue_names[link_slice],base.residue_names))
			atom_names = np.concatenate((link.atom_names[link_slice],base.atom_names))
			residue_indices = np.concatenate((link.residue_indices[link_slice],base.residue_indices))
			another = GMXStructure(pts=points,box=base.box,residue_indices=residue_indices,
				residue_names=residue_names,atom_names=atom_names)
			another.write(state.here+'monomer_%d_cg.gro'%(mono_num+1))
			backname = 'monomer_%d_back_aa'%(mono_num+1)
			#---backmap the combined monomer to make the topology using martinize
			#---! obviously this is a terrible hack
			#---! inconsistent suffixes below
			backmapper_wrapper('monomer_%d_cg.gro'%(mono_num+1),backname,cwd=state.here,jitter=True)
			#---make the topology
			gmx_run(state.gmxpaths['editconf']+' -f %s.gro -o %s.pdb'%(backname,backname),
				log='editconf-convert-pdb')
			monomer_name = 'monomer_%d'%(mono_num+1)
			topnames = self.martinize(name=monomer_name,structure=backname,
				#---override the final secondary structure
				secondary_structure=settings.complete_secondary_structure)
			#---! no redundancy check here is this inefficient?
			#---! note the following names are hardcoded assuming a single chain and hence Protein.itp
			#---! ...is output by martinize and then parsed by our wrapper around martinize 
			#---! ...so that the the names are e.g. monomer_1_protein
			if len(topnames)!=1 or len(self.contents[topnames[0]].molecules.keys())!=1:
				raise Exception('development. we found redundant topologies coming from martinize')
			else: 
				#---rename the molecule
				topname = topnames[0]
				molname = self.contents[topname].molecules.keys()[0]
				molname_new = 'monomer_%d'%(mono_num+1)
				self.contents[topname].molecules = {molname_new:
					self.contents[topname].molecules[molname]}
				self.contents[topname].molecules[molname_new]['moleculetype']['molname'] = molname_new
				new_itp = 'monomer_%d.itp'%(mono_num+1)
				self.contents[topname].write(state.here+new_itp)
				itp_names_equilibrate.append(new_itp)
				#---write the bonded version
				new_itp_bonded = 'monomer_%d_bonded.itp'%(mono_num+1)
				self.contents[topname].write(state.here+new_itp_bonded)
				#---we save the bonded version to the itp_names for now, but we will switch to the 
				#---...constraints version after minimization
				itp_names.append(new_itp_bonded)
				state.protein_prepared['composition'].append((molname_new,1))
		#---update the topology, mimicking the end of the original martinize wrapper at martini/martini.py
		state.itp = state.martinize_itps = itp_names
		state.martinize_itps_equilibrate = itp_names_equilibrate
		state.composition = state.protein_prepared['composition']
		vacuum_struct1 = GMXStructure(state.here+'monomer_1.gro')
		#---! temporary
		if True:
			vacuum_struct2 = GMXStructure(state.here+'monomer_2.gro')
			vacuum_struct1.add(vacuum_struct2)
		else: state.composition = [('monomer_1',1)]
		#---! center the protein which is important for later !!!
		vacuum_struct1.points -= vacuum_struct1.points.mean(axis=0)
		vacuum_struct1.write(state.here+'vacuum-unsized.gro')
		gmx('editconf',structure='vacuum-unsized',gro='vacuum',d=2,c=True,log='editconf-vacuum-box')
		write_topology('vacuum.top')
