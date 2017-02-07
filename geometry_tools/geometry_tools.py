#!/usr/bin/env python

"""
GEOMETRY TOOLS FOR AUTOMACS

Generic functions for basic three-dimensional construction tasks in various automacs codes.
"""

import numpy as np

_not_reported = ['principal_axis','plane_project','vecangle']

def principal_axis(pts):
	"""
	Return one of the three principal axes of a three-dimensional collection of points.
	"""
	axis = vecnorm(pts[0]-pts[-1])
	eigs = np.linalg.eig(np.dot(pts.T,pts))
	principal_axis_index = np.argsort(eigs[0])[-1]
	axis = vecnorm(eigs[1][:,principal_axis_index])
	return axis

def plane_project(x,n): 
	"""
	Project a vector x onto a plane normal to a vector n.
	"""
	return x-np.dot(x,n)/np.linalg.norm(n)*vecnorm(n)

def vecangle(v1,v2,degrees=False):
	"""
	Compute the angle between two vectors
	"""
	v1n,v2n = vecnorm(v1),vecnorm(v2)
	dotp = np.dot(v1n,v2n)
	angle = np.arccos(dotp)*(180./np.pi)
	if np.isnan(angle): return (0.0 if (v1n==v2n).all() else np.pi*(180/np.pi))
	if degrees: return angle
	else: return angle/180.*np.pi
