{

'table':{
#####
####
###
##
#
#---mimic a user copy command for realistic test sets
'tags':['aamd_cgmd','tag_support','tested_2017.09.18'],
'quick':"""
import amx,shutil
shutil.copyfile(amx.settings.ready,amx.settings.store)
""",
},

'bilayer_288_demo':{
#####
####
###
##
#
#---a legacy test used to make a small bilayer for the structure-repo
'tags':['cgmd','tested_2017.09.14'],
'metarun':[
{'step':'bilayer','do':'bilayer_control_cgmd','settings':"""
#---this demo generates a small coarse-grained bilayer for use in a test set
#---this test was used to generate @structure-repo/bilayers-cgmd/bilayer-cgmd-288.gro
step: bilayer
monolayer top: 144
composition top: {'DOPC':0.64,'DOPS':0.16,'POP2':0.2}
composition bottom: {'POPC':1.0}
thickness: 18
"""},
{'quick':'table','settings':"""
ready: s01-bilayer/md.part0001.gro
store: inputs/bilayer-cgmd-288.gro
"""},
]},

###---DEVELOPMENT NOTES
'comment_extras_testset':{'comment':"""

TESTSETS:
1. testset_bilayer_protein_free: attach helix0 structure from @structure-repo to a new, free bilayer
	requires 9.4min and is generally stable thanks to npt-bilayer equilibration bilayer_protein_adhesion
	may be less stable on dramatically different system sizes (recommend multiply step for large systems)
2. testset_bilayer_protein_flat: attach helix0 structure from @structure-repo to a new, flat bilayer
	equivalent to testset_bilayer_protein_free with added restraints
	note that users should use script-continue.sh to run the simulation until satisfactory binding
	restraints can be released with "make go bilayer_release" which was tested with the enth_demo experiment
3. testset_ultra1: a combination testset that includes items above. useful only for validating automacs

NOTES:
-- the testsets are somewhat slower than other examples (e.g. enth_demo) because they make new bilayers
-- the "table" step simulates a user who made one simulation and copied the result to inputs for another
-- no test sets to date (2017.09.18) work without a pre-made, *complete* protein structure
-- users who adapt these methods should be careful to check their topology and protein placement

"""},

'testset_bilayer_protein_free':{
#####
####
###
##
#
'tags':['cgmd','tested_2017.09.20','note_structure_repo_protein'],
'metarun':[
{'step':'bilayer','do':'bilayer_control_cgmd','settings':"""
step: bilayer
monolayer top: 90
monolayer bottom: 90
composition top: {'DOPC':0.64,'DOPS':0.16,'POP2':0.2}
composition bottom: {'POPC':1.0}
"""},
{'quick':'table','settings':"""
ready: s01-bilayer/md.part0001.gro
store: inputs/bilayer-cgmd-small.gro
"""},
{'step':'protein','do':'martinize','settings':"""
start structure: inputs/structure-repo/proteins/helix0.pdb
"""},
{'step':'adhere','do':'bilayer_protein_adhesion','settings':"""
force field: martini_upright_alt
sources: ['@martini/auto_ff/martini_upright_alt.ff']
placement method: banana
group up: resid 19
group down: resid 7
group origin: resid 7
bilayer structure: inputs/bilayer-cgmd-small.gro
protein_lattice:|{
	'nrows':1,'ncols':1,
	'lattice_type':'square',
	'space_scale':20,
	'total_proteins':1,
	'protein_shift_up':1.0,}
"""},
]},

'testset_bilayer_protein_flat':{
#####
####
###
##
#
'tags':['cgmd','tested_2017.09.20','note_structure_repo_protein'],
'metarun':[
{'step':'bilayer','do':'bilayer_control_flat','settings':"""
step: bilayer
monolayer top: 90
monolayer bottom: 90
composition top: {'DOPC':0.64,'DOPS':0.16,'POP2':0.2}
composition bottom: {'POPC':1.0}
"""},
{'quick':'table','settings':"""
ready: s01-bilayer/md.part0001.gro
store: inputs/bilayer-cgmd-small.gro
"""},
{'step':'protein','do':'martinize','settings':"""
start structure: inputs/structure-repo/proteins/helix0.pdb
"""},
{'step':'adhere','do':'bilayer_protein_adhesion','settings':"""
force field: martini_upright_alt
sources: ['@martini/auto_ff/martini_upright_alt.ff']
placement method: banana
group up: resid 19
group down: resid 7
group origin: resid 7
bilayer structure: inputs/bilayer-cgmd-small.gro
protein_lattice:|{
	'nrows':1,'ncols':1,
	'lattice_type':'square',
	'space_scale':20,
	'total_proteins':1,
	'protein_shift_up':1.0,}

#---EQUILIBRATION
equilibration: ['npt-bilayer-short','npt-bilayer']
mdp specs:|{
    'group':'cgmd',
    'mdps':{
        'input-em-steep-in.mdp':['minimize'],
        'input-md-npt-bilayer-short-eq-in.mdp':[{'restrain':'posre-com-only','pressure':'npt-semiisotropic-weak',
            'nsteps':500000,'groups':'protein','temperature':'protein','dt':0.001}],
        'input-md-npt-bilayer-eq-in.mdp':[{'restrain':'posre-com-only','pressure':'npt-semiisotropic-weak',
            'nsteps':500000,'groups':'protein','temperature':'protein','dt':0.01}],
        'input-md-in.mdp':[{'restrain':'posre-com-only','pressure':'npt-semiisotropic-weak',
            'nsteps':500000,'groups':'protein','temperature':'protein'}],},}

"""},
]},

'testset_ultra1':{
#####
####
###
##
#
#---a combination testset
'tags':['cgmd','dev'],
'metarun':[
{'step':'bilayer','do':'bilayer_control_cgmd','settings':"""
step: bilayer
monolayer top: 90
composition top: {'DOPC':0.64,'DOPS':0.16,'POP2':0.2}
composition bottom: {'POPC':1.0}
"""},
{'quick':'table','settings':"""
ready: s01-bilayer/md.part0001.gro
store: inputs/bilayer-cgmd-small.gro
"""},
{'step':'bilayer','do':'bilayer_control_flat','settings':"""
step: bilayer
monolayer top: 90
composition top: {'DOPC':0.64,'DOPS':0.16,'POP2':0.2}
composition bottom: {'POPC':1.0}
"""},
{'quick':'table','settings':"""
ready: s01-bilayer/md.part0001.gro
store: inputs/bilayer-cgmd-small-flat.gro
"""},
{'step':'protein','do':'martinize','settings':"""
start structure: inputs/helix0.pdb
"""},
{'step':'adhere','do':'bilayer_protein_adhesion','settings':"""
force field: martini-sources
sources: ['@martini/martini-sources.ff']
placement method: banana
group up: resid 19
group down: resid 7
group origin: resid 7
bilayer structure: inputs/bilayer-cgmd-small.gro
protein_lattice:|{
	'nrows':1,'ncols':1,
	'lattice_type':'square',
	'space_scale':20,
	'total_proteins':1,
	'protein_shift_up':1.0,}
"""},
{'step':'protein','do':'martinize','settings':"""
start structure: inputs/helix0.pdb
"""},
{'step':'adhere','do':'bilayer_protein_adhesion','settings':"""
force field: martini_upright_alt
sources: ['@martini/auto_ff/martini_upright_alt.ff']
placement method: banana
group up: resid 19
group down: resid 7
group origin: resid 7
bilayer structure: inputs/bilayer-cgmd-small-flat.gro
protein_lattice:|{
	'nrows':1,'ncols':1,
	'lattice_type':'square',
	'space_scale':20,
	'total_proteins':1,
	'protein_shift_up':1.0,
	}
"""},
]},

}