{

'helixer':{
#####
####
###
##
#
'tags':['cgmd'],
'script':'scripts/coiled_coil_linker_dev.py',
'params':'@bilayers/parameters.py',
'extensions':['geometry_tools/*.py','@bilayers/codes/bilayer.py'],
'settings':"""

USAGE NOTE:|
	developing helix routines
	...

step: helix
sequence:|
	MIPPQEASARRREIEDKLKQEEETLSFIRDSLEKSDQLTKNMVSILSSFESRLMKLENSI
	IPVHKQTENLQRLQENVEKTLSCLDHVISYYHVASDTEKIIREGPTGRLEEYLGSMAKIQ
	KAVEYFQDNSPDSPELNKVKLLFERGKESLESEFRSLMTRHSKVISPVLVLDLISADDEL
	EVQEDVVLEHLPESVLQDVIRISRWLVEYGRNQDFMNVYYQIRSSQLDRSIKGLKEHFRK
	SSSSSGVPYSPAIPNKRKDTPTKKPIKRPGRDDMLDVETDAYIHCVSAFVRLAQSEYQLL
	MGIIPEHHQKKTFDSLIQDALDGLMLEGENIVSAARKAIIRHDFSTVLTVFPILRHLKQT
	KPEFDQVLQGTAASTKNKLPGLITSMETIGAKALEDFADNIKNDPDKEYNMPKDGTVHEL
	TSNAILFLQQLLDFQETAGAMLASQETSSSATSYNSEFSKRLLSTYICKVLGNLQLNLLS
	KSKVYEDPALSAIFLHNNYNYILKSLEKSELIQLVAVTQKTAERSYREHIEQQIQTYQRS
	WLKVTDYIAEKNLPVFQPGVKLRDKERQMIKERFKGFNDGLEELCKIQKAWAIPDTEQRD
	KIRQAQKSIVKETYGAFLHRYSSVPFTKNPEKYIKYRVEQVGDMIDRLFDTSA
secondary structure: ''.join(['H' for i in range(1,90+1)])

#---assume an ideal coiled coil
oligomer specs:{'n_coils':2,'anti':True,'residue_shift':0}

residue_library: @martini/library-residues
sources: ['@martini/martini-sources.ff']
files: ['@martini/library-general-structs/martini-water.gro']
force field: martini-sources
sequence slice: "1-90"
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

'helixer_old_settings':{
#####
####
###
##
#
'tags':[],
'script':'scripts/structure_explorer.py',
'params':'parameters.py',
'extensions':['geometry_tools/structural_biology.py'],
'settings':"""

USAGE NOTE:|
	developing helix routines
	...

step: helix
sequence: |
	MIPPQEASARRREIEDKLKQEEETLSFIRDSLEKSDQLTKNMVSILSSFESRLMKLENSI
	IPVHKQTENLQRLQENVEKTLSCLDHVISYYHVASDTEKIIREGPTGRLEEYLGSMAKIQ
	KAVEYFQDNSPDSPELNKVKLLFERGKESLESEFRSLMTRHSKVISPVLVLDLISADDEL
	EVQEDVVLEHLPESVLQDVIRISRWLVEYGRNQDFMNVYYQIRSSQLDRSIKGLKEHFRK
	SSSSSGVPYSPAIPNKRKDTPTKKPIKRPGRDDMLDVETDAYIHCVSAFVRLAQSEYQLL
	MGIIPEHHQKKTFDSLIQDALDGLMLEGENIVSAARKAIIRHDFSTVLTVFPILRHLKQT
	KPEFDQVLQGTAASTKNKLPGLITSMETIGAKALEDFADNIKNDPDKEYNMPKDGTVHEL
	TSNAILFLQQLLDFQETAGAMLASQETSSSATSYNSEFSKRLLSTYICKVLGNLQLNLLS
	KSKVYEDPALSAIFLHNNYNYILKSLEKSELIQLVAVTQKTAERSYREHIEQQIQTYQRS
	WLKVTDYIAEKNLPVFQPGVKLRDKERQMIKERFKGFNDGLEELCKIQKAWAIPDTEQRD
	KIRQAQKSIVKETYGAFLHRYSSVPFTKNPEKYIKYRVEQVGDMIDRLFDTSA
secondary structure: ''.join(['H' for i in range(1,90+1)])
sources: ['@martini/martini-sources.ff']
spacing: 0.35
sequence slice: "1-90"
martinize path: '~/libs/martinize.py'
residue_library: 'inputs/martini/library-residues'
backmapper path: '~/libs/martini-backmapper'
crude structure: helix-crude
atomistic structure: helix
atomistic minimization: true
atomistic force field: charmm36
force field: martini
ff includes:| ['martini-v2.2','martini-v2.0-lipids',
 'martini-v2.2-aminoacids','martini-v2.0-ions']
water: tip3p
mdp specs atomistic:|
 {'group':'aamd',
 'input-em-steep-in.mdp':['minimize',{'nsteps':10000}],
 'input-md-in.mdp':None}
mdp specs:|
 {'group':'cgmd',
 'input-em-steep-in.mdp':['minimize'],
 'input-md-in.mdp':[{'temperature':'protein-vacuum','dt':0.001,
 'groups':'blank','pressure':'npt-isotropic-weak','nsteps':100000,'nstxtcout':1000}],
 'input-md-npt-eq-in.mdp':[{'temperature':'protein-vacuum','dt':0.001,
 'groups':'blank','pressure':'npt-isotropic-weak','nsteps':100000,'nstxtcout':1000}]}
equilibration: 'npt'
spacecurve: 'lambda x,pitch=1.0,yaw=0.1:(pitch*np.cos(yaw*x),pitch*np.sin(yaw*x),x)'
lattice: [[0,0,0],[0,0,2.5]]
box buffer: 10.0
water box: @martini/library-general-structs/martini-water.gro
files: []
atom resolution: cgmd
sol: W
run part two: True
orientations: [0,180]
wind_backwards: False
ionic strength: 0.150
cation: NA+
anion: CL-
protein water gap: 3

"""},

}