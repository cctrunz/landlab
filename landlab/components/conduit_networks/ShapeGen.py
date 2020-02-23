from numpy import sin, cos, pi, sqrt, fabs, arctan2, linspace
import numpy as np


# Generate x, y points for an ellipse. Can set rotation angle with theta
def genEll(r1,r2,theta=0,d=1000):

	t=linspace(0, 2*pi-2*pi/d, d-1)

	#If we pass arrays of r1 and r2, then vectorize the calculation of xc shapes
	if len(r1)>1:
		r1 = np.broadcast_to(r1, (len(t),len(r1)) ).transpose()
		r2 = np.broadcast_to(r2, (len(t),len(r2)) ).transpose()
		t = np.broadcast_to(t, (len(r1),len(t)) )

	x=r1*cos(t)
	y=r2*sin(t)


	if(theta!=0):

		tx=x
		ty=y

		x=tx*cos(theta)-ty*sin(theta)
		y=tx*sin(theta)+ty*cos(theta)

	return x, y

# Generate x, y points for a circle
def genCirc(r,d=1000):

	return genEll(r, r,d=d)

# Generate x, y points for the lower half of a circle
def genSemiCirc(r, d=1000):

	return genSemiEll(r, r,d=d)

# Generate x, y points for the lower half of an ellipse
def genSemiEll(r1,r2,d=1000):

	t=linspace(-1*pi, 0, d)
	x=r1*cos(t)
	y=r2*sin(t)

	return x, y
