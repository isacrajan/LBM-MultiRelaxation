#
"""
Start Date: 14 August 2017
@author: Isac Rajan
Single Sided Lid Driven Cavity using LBM MRT Scheme
Email:	isacrajan@gmail.com
"""
import numpy as np

#No of grids in X and Y direction
NX = 300
NY = 300

#Defining the parameters
TSTEPS = 150000
cx = np.array([0,1,0,-1,0,1,-1,-1,1])
cy = np.array([0,0,1,0,-1,1,1,-1,-1])
w = np.array([4./9,1./9,1./9,1./9,1./9,1./36,1./36,1./36,1./36])
uo = 0.1 #Lid Velocity
rhoo = 0.1
Re = 100 #Reynolds Number
nu = uo*NY/Re #Kinematic Viscosity
omega = 1.0/(3.0*nu + 0.5) #inverse of tau, as in formulation

#Initializing the arrays
u = np.zeros((NX,NY))
v = np.zeros((NX,NY))
rho = rhoo * np.ones((NX,NY))
f = np.zeros((9,NX,NY))
feq = np.zeros((9,NX,NY))
m = np.empty((9,NX,NY))
meq = np.empty((9,NX,NY))
M = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1],
							[-4, -1, -1, -1, -1, 2, 2, 2, 2],
							[4, -2, -2, -2, -2, 1, 1, 1, 1],
							[0, 1, 0, -1, 0, 1, -1, -1, 1],
							[0, -2, 0, 2, 0, 1, -1, -1, 1],
							[0, 0, 1, 0, -1, 1, 1, -1, -1],
							[0, 0, -2, 0, 2, 1, 1, -1, -1],
							[0, 1, -1, 1, -1, 0, 0, 0, 0],
							[0, 0, 0, 0, 0, 1, -1, 1, -1]]) #Transformation Matrix
Minv = np.linalg.inv(M) #Inverse of Transformation Matrix							
S = np.diag([0,omega,omega,0,omega,0,omega,omega,omega]) #Relaxation Matrix

#Initialization of Lid Velocity -- Top Lid
u[:][NY-1] = uo;
v[:][NY-1] = 0.0;

#Evaluating feq at t=0 & Initializing f
t1 = np.multiply(u,v)
for k in range(9):
	t2 = np.add(np.multiply(u,cx[k]), np.multiply(v, cy[k]))
	feq[k][:][:] = rho[:][:] * w[k] * (1.0 + 3.0*t2 + 4.5*t2*t2 - 1.50*t1)
	f[k][:][:] = feq[k][:][:]

#Evaluating meq at t=0 & Initializing m
for k in range(9):
	if k == 0:
		meq[k][:][:] = rho[:][:]
	elif k == 1:
		meq[k][:][:] = -2*rho[:][:] + 3*(np.add(np.multiply(u,u),np.multiply(v,v)))
	elif k == 2:
		meq[k][:][:] = rho[:][:] - 3*(np.add(np.multiply(u,u),np.multiply(v,v)))
	elif k == 3:
		meq[k][:][:] = np.multiply(rho,u)
	elif k == 4:
		meq[k][:][:] = -u[:][:]
	elif k == 5:
		meq[k][:][:] = np.multiply(rho,v)
	elif k == 6:
		meq[k][:][:] = -v[:][:]
	elif k == 7:
		meq[k][:][:] = np.multiply(u,u) - np.multiply(v,v) 
	elif k == 8:
		meq[k][:][:] = 	np.multiply(u,v)
	m[k][:][:] = meq[k][:][:]







