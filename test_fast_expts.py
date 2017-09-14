{

'table':{
#####
####
###
##
#
'tags':['aamd_cgmd','tag_support'],
'quick':"""
import amx,shutil
shutil.copyfile(amx.settings.ready,amx.settings.store)
""",
},

'ultra1':{
#####
####
###
##
#
'tags':['cgmd','tag_?!testset'],
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
bilayer structure: inputs/structure-repo/bilayers-cgmd/bilayer-cgmd-small-flat.gro
protein_lattice:|{
	'nrows':1,'ncols':1,
	'lattice_type':'square',
	'space_scale':20,
	'total_proteins':1,
	'protein_shift_up':1.0,
	}
"""},
]},

'bilayer_demo_288':{
#####
####
###
##
#
'tags':['cgmd','tested_2017.09.14','tag_continues'],
'metarun':[
{'step':'bilayer','do':'bilayer_control_cgmd','settings':"""
step: bilayer
monolayer top: 144
composition top: {'DOPC':0.64,'DOPS':0.16,'POP2':0.2}
composition bottom: {'POPC':1.0}
thickness: 18
"""},
{'quick':'table','settings':"""
ready: s01-bilayer/md.part0001.gro
store: inputs/structure-repo/bilayers-cgmd/bilayer-cgmd-288.gro
"""},
]},

'enth_demo':{
#####
####
###
##
#
'tags':['cgmd','tag_structure_repo'],
'metarun':[
{'step':'protein','do':'martinize','settings':"""
start structure: inputs/structure-repo/proteins/1H0A-prepped.pdb
"""},
{'step':'adhere','do':'bilayer_protein_adhesion','settings':"""
force field: martini_upright_alt
sources: ['@martini/auto_ff/martini_upright_alt.ff']
placement method: globular_up_down
group up: all
group down: ['resid 1-22','resid 68-72']
group origin: ['resid 1-22','resid 68-72']
bilayer structure: inputs/structure-repo/bilayers-cgmd/bilayer-cgmd-288.gro
protein_lattice:|{
	'nrows':1,'ncols':1,
	'lattice_type':'square',
	'space_scale':20,
	'total_proteins':1,
	'protein_shift_up':1.0,
	}
"""},
]},

'ultra2':{
#####
####
###
##
#
'tags':['cgmd','tag_?!testset'],
'metarun':[
{'step':'bilayer','do':'bilayer_control','settings':"""
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
	'protein_shift_up':1.0,
	}
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