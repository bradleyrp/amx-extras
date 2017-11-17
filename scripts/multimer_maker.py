#!/usr/bin/env python

"""
MULTIMER MAKER
one-stop shop for complicated proteins
"""

from amx import *
init()
make_step(settings.step)
write_mdp()
ProteinMultimer()
minimize('vacuum')
solvate(
	structure='vacuum-minimized',
	gro='solvate')
write_top('solvate.top')
minimize('solvate')
#---extract the protein from the water
gmx('make_ndx',structure='vacuum',
	#---! hack below to get the protein out
	inpipe='keep 1\nq\n',ndx='index-protein',log='group-protein')
#---extract the protein
gmx('trjconv',structure='solvate-minimized',gro='protein-ready-uncentered',tpr='em-solvate-steep',
	ndx='index-protein',log='trjconv-extract-protein')
#---center the protein at the origin for the adhesion routine
gmx('editconf',structure='protein-ready-uncentered',gro='protein-ready',
	center="0 0 0",log='editconf-recenter')
#---prepare for the subsequent adhesion step
state.protein_prepared = {'gro':state.here+'protein-ready.gro',
	#---! clumsy handling of ITP here and additions to adhere_protein_bilayer.py
	state.here+'top':'vacuum.top','itps':[os.path.join(state.here,i) for i in state.itp],
	#---! currently hardcoded
	'protein_composition':[('monomer_1',1),('monomer_2',1)]}
#---the bilayer adhesion either takes one ITP for atomistic or martinize_itps for MARTINI
#---! overall that needs to be more carefully reworked
state.martinize_itps = state.itps
