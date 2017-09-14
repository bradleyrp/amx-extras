{

'multiply':{
#####
#---note that this was tested with e.g. bilayer_control_multiply
'tags':['aamd_cgmd','tested_2017.09.13'],
'script':'scripts/multiply.py',
'params':None,
'extensions':['*.py','../bilayers/codes/*.py'],
'settings':"""
step: large
requires: multiply
equilibration: npt-bilayer
proceed: True
genconf gap: 0.3
nx: 10
ny: 10
"""
},

}