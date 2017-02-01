#!/usr/bin/env python

"""
MESH ROUTINES

Useful for making triangulated surfaces.
"""

import numpy as np
import scipy
import scipy.spatial
import scipy.interpolate
from numpy import linalg

_not_reported = ['triarea','vecnorm','facenorm','torusnorm','reclock','beyonder','makemesh','rotation_matrix']
_shared_extensions = ['vecnorm','rotation_matrix','makemesh']

#---geometry functions
triarea = lambda a : linalg.norm(np.cross(a[1]-a[0],a[2]-a[0]))/2.
vecnorm = lambda vec: np.array(vec)/linalg.norm(vec)
facenorm = lambda a: np.cross(a[1]-a[0],a[2]-a[0])

def rotation_matrix(axis,theta):
	"""
	Return the rotation matrix associated with counterclockwise rotation about
	the given axis by theta radians using Euler-Rodrigues formula.
	"""
	axis = np.asarray(axis)
	theta = np.asarray(theta)
	if all(axis==0): return np.identity(3) 
	axis = axis/np.sqrt(np.dot(axis,axis))
	a = np.cos(theta/2)
	b, c, d = -axis*np.sin(theta/2)
	aa, bb, cc, dd = a*a, b*b, c*c, d*d
	bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
	return np.array([[aa+bb-cc-dd,2*(bc+ad),2*(bd-ac)],[2*(bc-ad),aa+cc-bb-dd,2*(cd+ab)],
		[2*(bd+ac),2*(cd-ab),aa+dd-bb-cc]])

def reclock(points,ords=[0,1,2]):
	"""
	Reorders R2 points in clockwise order.
	"""
	rels = points[ords]-np.mean(points[ords],axis=0)
	return np.argsort(np.arctan2(rels[:,0],rels[:,1]))[::-1]

def beyonder(points,vecs,dims=[0,1],growsize=0.2,new_only=False,growsize_nm=None,return_ids=False):
	"""
	Given points and box vectors, this function will expand a bilayer patch across periodic boundaries. It is 
	meant to create inputs for meshing algorithms which later store data under the correct topology. The dims
	argument names the wrapping dimensions, which is typically XY.
	"""
	over = np.array([(i,j,k) for i in [-1,0,1] for j in [-1,0,1] for k in [-1,0,1]])
	over = np.array([list(j) for j in list(set([tuple(i) for i in over[:,:len(dims)]])) 
		if j != tuple(np.zeros(len(dims)))])
	vec2 = np.array([(vecs[i] if i in dims else 0.) for i in range(3)])	
	over2 = np.array([[(o[i] if i in dims else 0.) for i in range(3)] for o in over])
	alls = np.concatenate([points]+[points+i*vec2 for i in over2])
	#---save indices for extra points
	inds = np.concatenate([np.arange(len(points)) for i in range(len(over2)+1)])
	#---bugfix applied 2015.07.20 to fix incorrect ghost index numbering
	if 0: inds = np.concatenate([np.arange(len(points))]+[np.ones(len(over))*j 
		for j in np.arange(len(points))])
	#---override the proportional growsize by a literal size
	if growsize_nm != None: growsize = max([growsize_nm/vecs[i] for i in dims])
	valids = np.where(np.all([np.all((alls[:,d]>-1*vec2[d]*growsize,alls[:,d]<vec2[d]*(1.+growsize)),axis=0) 
		for d in dims],axis=0))[0]
	if not return_ids: return alls[valids]
	else: return alls[valids],inds[valids].astype(int)

def torusnorm(pts1,pts2,vecs):

	"""
	Compute distances between points on a torus.
	"""

	cd = np.array([scipy.spatial.distance.cdist(pts1[:,d:d+1],pts2[:,d:d+1]) for d in range(2)])
	cd[0] -= (cd[0]>vecs[0]/2.)*vecs[0]
	cd[1] -= (cd[1]>vecs[1]/2.)*vecs[1]
	cd2 = linalg.norm(cd,axis=0)
	return cd2

def makemesh(pts,vec,growsize=0.2,curvilinear_neighbors=10,
	curvilinear=True,debug=False,growsize_nm=None,excise=True,areas_only=False):

	"""
	Function which computes curvature and simplex areas on a standard mesh.
	"""

	if debug: 
		import time
		st = time.time()
		def checkpoint(note):
			print(note)
			st = time.time()

	nmol = len(pts)
	pts = pts
	vec = vec
	if debug: 
		import time
		st = time.time()
		print("[STATUS] start makemesh %0.2f"%(time.time()-st))
	ptsb,ptsb_inds = beyonder(pts,vec,growsize=growsize,growsize_nm=growsize_nm,return_ids=True)
	if debug: print("[STATUS] project curvilinear="+str(curvilinear)+" %0.2f"%(time.time()-st))
	#---if curvilinear then use the isomap otherwise project onto the xy plane
	if curvilinear: 
		import sklearn
		from sklearn import manifold
		proj = manifold.Isomap(curvilinear_neighbors,2).fit_transform(ptsb)
	else: proj = ptsb[...,:2]
	if debug: checkpoint("[STATUS] delaunay %0.2f"%(time.time()-st))
	if debug: checkpoint("[STATUS] shape="+str(np.shape(ptsb)))
	dl = scipy.spatial.Delaunay(proj)
	if debug: checkpoint("[STATUS] reclock %0.2f"%(time.time()-st))
	simplices = np.array([a[reclock(ptsb[a])] for a in dl.simplices])
	#---rework simplices and ptsb to exclude superfluous points
	if debug: checkpoint("[STATUS] trim %0.2f"%(time.time()-st))
	#---relevants is a unique list of simplices with exactly one member that is equal to a core vertex point
	relevants = np.unique(np.concatenate([simplices[
		np.where(np.sum(simplices==i,axis=1)==1)[0]] for i in range(nmol)]))
	points = ptsb[relevants]
	ghost_indices = ptsb_inds[relevants]
	ptsb = points
	if debug: checkpoint("[STATUS] simplices %0.2f"%(time.time()-st))
	simplices = np.array([[np.where(relevants==r)[0][0] for r in s] 
		for s in simplices if np.all([r in relevants for r in s])])
	#---end rework
	if debug: checkpoint("[STATUS] areas %0.2f"%(time.time()-st))
	areas = np.array([triarea(ptsb[a]) for a in simplices])
	if areas_only: return {'simplices':simplices,'areas':areas,'nmol':nmol,'vec':vec,'points':points}
	if debug: checkpoint("[STATUS] facenorms %0.2f"%(time.time()-st))
	facenorms = np.array([vecnorm(facenorm(ptsb[a])) for a in simplices])	
	if debug: checkpoint("[STATUS] vertex-to-simplex %0.2f"%(time.time()-st))
	v2s = [np.where(np.any(simplices==i,axis=1))[0] for i in range(nmol)]
	if debug: checkpoint("[STATUS] vertex normals %0.2f"%(time.time()-st))
	vertnorms = np.array([vecnorm(np.sum(facenorms[ind]*\
		np.transpose([areas[ind]/np.sum(areas[ind])]),axis=0)) for ind in v2s])
	principals = np.zeros((nmol,2))
	nl = []
	if debug: checkpoint("[STATUS] curvatures %0.2f"%(time.time()-st))
	for v in range(nmol):
		neighbors = np.unique(simplices[np.where(np.any(simplices==v,axis=1))[0]])
		neighbors = neighbors[neighbors!=v]
		nl.append(neighbors)
		edges = ptsb[neighbors]-ptsb[v]
		weights = [areas[sl]/2./np.sum(areas[v2s[v]]) for sl in v2s[v]]
		tijs = [vecnorm(np.dot(np.identity(3)-np.outer(vertnorms[v],
			vertnorms[v].T),ab)) for ab in edges]
		kijs = [np.dot(vertnorms[v].T,ab)/linalg.norm(ab)**2 for ab in edges]
		ct = np.sum([weights[ind]*kijs[ind]*np.outer(tijs[ind],tijs[ind]) 
			for ind,i in enumerate(v2s[v])],axis=0)
		wsign = 1-2*(linalg.norm(np.array([1,0,0])+\
			vertnorms[v])<linalg.norm(np.array([1,0,0])-vertnorms[v]))
		wvi = vecnorm(np.array([1,0,0])+wsign*vertnorms[v])
		hm = np.identity(3)-2*np.outer(wvi,wvi.T)
		hhm = np.dot(np.dot(hm.T,ct),hm)
		principals[v] = -1*hhm[1,1],-1*hhm[2,2]
	if debug: checkpoint("[STATUS] PBC neighborlist %0.2f"%(time.time()-st))
	#---neighborlist under PBCs
	checksubssort,nlsubs = np.where(torusnorm(points[nmol:],points[:nmol],vec)==0)
	#if not all(checksubssort==np.arange(len(points)-nmol)): raise Exception('torusnorm lookup fail')
	try: nlpbc = [[(i if i<nmol else nlsubs[i-nmol]) for i in n] for n in nl]
	except: nlpbc = []
	gauss = (3*principals[:,0]-principals[:,1])*(3*principals[:,1]-\
		principals[:,0])
	mean = 1./2*((3*principals[:,0]-principals[:,1])+\
		(3*principals[:,1]-principals[:,0]))
	if debug: checkpoint("[STATUS] complete %0.2f"%(time.time()-st))
	return {'nmol':nmol,'vec':vec,'simplices':simplices,'points':points,
		'areas':areas,'facenorms':facenorms,'vertnorms':vertnorms,'principals':principals,
		'ghost_ids':ghost_indices,'gauss':gauss,'mean':mean}
