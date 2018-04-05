{

'multiply':{
#####
####
###
##
#
'tags':['aamd_cgmd','tested_2017.09.18'],
'script':'scripts/multiply.py',
'params':None,
'extensions':['*.py','../bilayers/codes/*.py'],
'settings':"""
step: large
requires: multiply
equilibration: npt-bilayer
maxwarn: 1
proceed: True
genconf gap: 0.3
nx: 2
ny: 2
"""
},

'multiply_general':{
#####
####
###
##
#
'tags':['aamd_cgmd','tested_2017.09.18'],
'script':'scripts/multiply-general.py',
'params':None,
'extensions':['*.py','../bilayers/codes/*.py'],
'settings':"""
step: large
requires: multiply
equilibration: npt-bilayer
maxwarn: 1
proceed: True
genconf gap: 0.3
nx: 2
ny: 2
"""
},

'demo_homology_helix0':{
#####
####
###
##
#
'tags':['aamd','test'],
'metarun':[
{'step':'homology','do':'homology','settings':"""
start structure: @structure-repo/proteins/helix0.pdb
template chain : {'A':{'startres':1,'stopres':22}}
point mutation : {'A':['M10T']}
numbering: cryst
modeller script: @homology/modeller_script.py
refinement: fast
"""}]},

}
