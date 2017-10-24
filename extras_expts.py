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
proceed: True
genconf gap: 0.3
nx: 2
ny: 2
"""
},

'homology_demo':{
#####
####
###
##
#
'tags':['aamd','tag_TEST?'],
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
