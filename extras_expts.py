{

'multiply':{
#####
'tags':[],
'script':'scripts/multiply.py',
'params':'parameters.py',
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