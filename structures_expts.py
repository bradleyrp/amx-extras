{

'helixer':{
	#####
	####
	###
	##
	#
	'tags':['cgmd','tag_exo70'],
	'script':'scripts/coiled_coil_linker_dev.py',
	'params':'@bilayers/parameters.py',
	'extensions':['geometry_tools/*.py','@bilayers/codes/bilayer.py'],
	'settings':"""

	USAGE NOTE:|
		this creates a coiled-coil structure (the "helixer" name is clumsy)

	step: helix
	sequence:|
		TIPFAPLSQSKTSHGVSLSIVMEGDYRDSQYTTSYPIDKRHGGSFSADTRQAVTGSESHT
		LKDPETNPLVLNREFTIQEDDKDPIQLKVTPLKDNKSRYYIGSFEPRCLALLFVEYTIAC
		PILLLQSRAGDVEVRRPEGIESGALLHHRMLSIMMVLQLDEPKLFSLPEKDYIKTFGDDS
		LISEVRLKQVAELAEILIGKNQYYETIQSFLQYVASEKRKLKVLKMERFAVGSMIASNIF
		EPLSQKKRLEVASLYYIKSFLCYIKNPYDELRRKADLHTIESAMKSAILSVQDVSGHDRT
		GDKFKKQDINNNEAPLSLVKFLQRMFEEDGYPNKEVESLERALKSNILKYQCSKFAYPTG
		TKDLVNADDKHKEEDQQDEKHADREQSTDLLEKTFVAHMVVESELLRSELFSQPKQRREN
		RSSGHELRKLPGKVAQSLQAVSFMDEHSIQTYVEIISLTKNNMAILNPLNENFDTKLSES
		SIATKLELRIEKILMAAGTIISPKNASDARVPGEQGKIIVLDAEDDSISELFYYKQTVNF
		SSDYRKLRKSIHTWETQSQFRQEQSLINRELVVVLITYIRALGFAKLQENDQDIRNSDME
		VLYDPTKLQVWVPTKKEWIHLKSQMEIVISITNEKTERRAQPQLVKEHLRPQQ
	secondary structure: ''.join(['H' for i in range(1,90+1)])
	sequence slice: "1-90"

	#---assume an ideal coiled coil
	oligomer specs:{'n_coils':2,'anti':True,'residue_shift':0}

	residue_library: @martini/library-residues
	sources: ['@martini/martini-sources.ff']
	files: ['@martini/library-general-structs/martini-water.gro']
	force field: martini-sources
	review 3d: false
	backwards script: @martini/bin/backward/backward.py
	solvent: martini-water
	martinize script: @martini/bin/martinize.py
	atomistic force field: charmm27
	ionic_strength: 0.150
	water: tip3p
	sol: W
	cation: NA+
	anion: CL-
	equilibration: short
	water edges: [2,2,5]

	mdp specs:|{
	    'group':'cgmd',
	    'mdps':{
	        'input-em-steep-in.mdp':['minimize'],
	        'input-md-short-eq-in.mdp':[{'nsteps':1000000,'groups':'protein-water',
	        	'temperature':'protein-water','dt':0.001}],
	        'input-md-in.mdp':[{'nsteps':10000000,'groups':'protein-water',
	        	'temperature':'protein-water','restrain':'posre-com-only','pressure':'standard-isotropic'}],
	        },
	    }

	"""},

'coiled_coil_only':{
	#####
	####
	###
	##
	#
	'tags':['cgmd','tag_exo70'],
	'script':'scripts/coiled_coil_only.py',
	'params':'@bilayers/parameters.py',
	'extensions':['geometry_tools/*.py','@bilayers/codes/bilayer.py'],
	'settings':"""

	USAGE NOTE:|
		forked from "helixer" to create the helix without solvent
		meant to precede a flat-bilayer routine

	step: helix
	sequence:|
		TIPFAPLSQSKTSHGVSLSIVMEGDYRDSQYTTSYPIDKRHGGSFSADTRQAVTGSESHT
		LKDPETNPLVLNREFTIQEDDKDPIQLKVTPLKDNKSRYYIGSFEPRCLALLFVEYTIAC
		PILLLQSRAGDVEVRRPEGIESGALLHHRMLSIMMVLQLDEPKLFSLPEKDYIKTFGDDS
		LISEVRLKQVAELAEILIGKNQYYETIQSFLQYVASEKRKLKVLKMERFAVGSMIASNIF
		EPLSQKKRLEVASLYYIKSFLCYIKNPYDELRRKADLHTIESAMKSAILSVQDVSGHDRT
		GDKFKKQDINNNEAPLSLVKFLQRMFEEDGYPNKEVESLERALKSNILKYQCSKFAYPTG
		TKDLVNADDKHKEEDQQDEKHADREQSTDLLEKTFVAHMVVESELLRSELFSQPKQRREN
		RSSGHELRKLPGKVAQSLQAVSFMDEHSIQTYVEIISLTKNNMAILNPLNENFDTKLSES
		SIATKLELRIEKILMAAGTIISPKNASDARVPGEQGKIIVLDAEDDSISELFYYKQTVNF
		SSDYRKLRKSIHTWETQSQFRQEQSLINRELVVVLITYIRALGFAKLQENDQDIRNSDME
		VLYDPTKLQVWVPTKKEWIHLKSQMEIVISITNEKTERRAQPQLVKEHLRPQQ
	secondary structure: ''.join(['H' for i in range(1,90+1)])
	sequence slice: "1-90"

	#---assume an ideal coiled coil
	oligomer specs:{'n_coils':2,'anti':True,'residue_shift':0}

	residue_library: @martini/library-residues
	sources: ['@martini/martini-sources.ff']
	files: ['@martini/library-general-structs/martini-water.gro']
	force field: martini-sources
	review 3d: false
	backwards script: @martini/bin/backward/backward.py
	martinize script: @martini/bin/martinize.py
	atomistic force field: charmm27

	mdp specs:|{
	    'group':'cgmd',
	    'mdps':{
	        'input-em-steep-in.mdp':['minimize'],
	        },
	    }

	"""},

}