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

}