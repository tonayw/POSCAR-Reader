import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import CubicSpline
import sys, getopt
import math

# Global variables
eps_0 = 8.8541878176E-12
q_e = 1.60217662E-19
A = 1E-10

# returns the conventional lattice information from primitive cell
def setPOSCARinfo(POSCAR, lattice, ionCHG):
	#read file and parse!
	file = open(POSCAR,"r")
	lines = file.read().split("\n")

	# set scaling factor
	sf = float(lines[1])

	# set primitive lattice vector
	pVec1 = list(map(float,lines[2].split()))
	pVec2 = list(map(float,lines[3].split()))
	pVec3 = list(map(float,lines[4].split()))

	# multiply scaling factor to the lattice vectors
	pVec1 = [e * sf for e in pVec1]
	pVec2 = [e * sf for e in pVec2]
	pVec3 = [e * sf for e in pVec3]

	# set lattice parameters
	a1 = np.linalg.norm(pVec1)
	a2 = np.linalg.norm(pVec2)
	a3 = np.linalg.norm(pVec3)

	# set atoms
	atoms = lines[5].split()
	atomN = list(map(int,lines[6].split()))
	N = 0
	for e in atomN:
		N = N + e

	# set atom ion charges
	ionchg = [[""] for e in range(N)]
	ii, jj = 0, 0
	for atom in atomN:
		for e in range(atom):
			ionchg[ii] = ionCHG[jj]
			# print("hahahahah")
			ii = ii + 1
		jj = jj + 1

	# set cell volume
	vol = np.dot(pVec3, np.cross(pVec1, pVec2))

	# set conventional lattice constant
	if(lattice[0] == "Z"):
		ratio = 2/np.sqrt(3)/np.sqrt(2) 	# a-to-h ratio
		lat_a, lat_c = a1*ratio*3**0.5, a1*ratio*3**0.5
		h = a3*ratio
	elif(lattice[0] == "W"):
		ratio = 1
		lat_a, lat_c = a1, a3
		h = a3
	elif(lattice[0] == "F"):
		ratio = 2/np.sqrt(3)/np.sqrt(2) 	# a-to-h ratio
		lat_a, lat_c = a1*ratio*3**0.5, a1*ratio*3**0.5
		h = a3*ratio
	elif(lattice[0] == "1"):
		ratio = math.sin(93.887*math.pi/180) 	# a-to-h ratio
		lat_a, lat_c = a1, a3
		h = a3*ratio
	else:
		ratio = 1
		lat_a, lat_c = a1, a3
		h = a3

	# print("lat_a, lat_c = %.4f, %.4f" % (lat_a, lat_c))
	# print("h = %.4f" % h)
	# print("vol = %.4f" % (vol))
	# print("ratio = %.4f" % (ratio))

	if(lines[7][0] == "D" or lines[7][0] == "d"):
		coord = "Direct"
	elif(lines[7][0] == "D" or lines[7][0] == "d"):
		coord = "Cart"
	else:
		print("POSCAR is neither Direct nor Cart. SOMETHING IS WRONG!!")
		coord = "SOMETHING WRONG"

	# set atomic coordinates
	crd = [[""] for e in range(N)]
	ii = 0
	for e in range(N):
		crd[ii] = list(map(float,lines[ii+8].split()[:3]))[2]*h
#		crd[ii] = list(map(float,lines[ii+8].split()))[2]*h
		ii = ii + 1

	# print(atomN)
	ion = [["",""] for e in range(N)]
	if(coord == "Direct"):
		for e, i, cr in zip(ion, ionchg, crd):
			e[0] = cr
			e[1] = i
	else:
		print("Cartesian coordinate not supported!!")
	return lat_a, lat_c, h, vol, ratio, ion

def getOriginalMultipoles(chg,ion,vol,ratio):
	# add one more data for periodic boundary condition
	chg.append([chg[-1][0] + chg[1][0], chg[0][1]])

	# set x and y
	x = np.asarray(list(zip(*chg))[0])
	y = np.asarray(list(zip(*chg))[1])

	# do this for charge nomalization regarding added data above
	y[0] = y[0]/2
	y[-1] = y[-1]/2

	# calculate multipoles with original data!
	M_i, M_e, M_e_new = 0, 0, 0
	D_i, D_e, D_e_new = 0, 0, 0
	Q_i, Q_e, Q_e_new = 0, 0, 0
	for e in ion:
		M_i = M_i + e[1]
		D_i = D_i + e[0]*A*e[1]
		Q_i = Q_i + (e[0]*A)**2*e[1]
	for xx, yy in zip(x, y):
		M_e = M_e - yy*x[1]
		D_e = D_e - xx*A*yy*x[1] 
		Q_e = Q_e - (xx*A)**2*yy*x[1]
	scaling = 1/(vol*A**3)*q_e*q_e/(2*eps_0)/q_e
	M_e  = M_e/ratio
	D_i = D_i*scaling*A
	D_e = D_e*scaling*A/ratio
	Q_i = Q_i*scaling
	Q_e = Q_e*scaling/ratio
	chg.pop(-1)
	return [[M_i, M_e, M_i + M_e],
			[D_i, D_e, D_i + D_e],
			[Q_i, Q_e, Q_i + Q_e]]

def normlizeChg(chg,ion,ratio):
	# set x and y
	x = np.asarray(list(zip(*chg))[0])
	y = np.asarray(list(zip(*chg))[1])

	# calculate multipoles with original data!
	tmp = 0
	M_i, M_e, M_e_new = 0, 0, 0
	for e in ion:
		M_i = M_i + e[1]
	for yy in y:
		M_e = M_e - yy*x[1]
	M_e  = M_e/ratio

	for yy in y:
		tmp = tmp - yy**2/ratio/ratio*x[1]
	for yy in y:
		M_e_new = M_e_new - yy*x[1] + yy**2/tmp*(M_i+M_e)*x[1]/ratio
	M_e_new = M_e_new/ratio
	# print(M_e_new)
	chgnorm = [["",""] for e in range(len(chg))]
	for e, en in zip(chg, chgnorm):
		en[0] = float(e[0])
		en[1] = float(e[1] - e[1]**2/tmp*(M_i+M_e)/ratio)

	return chgnorm

def getMonopole(chg,ion,grid,ratio):
	# add one more data for periodic boundary condition
	chg.append([chg[-1][0] + chg[1][0], chg[0][1]])

	# set x and y
	x = np.asarray(list(zip(*chg))[0])
	y = np.asarray(list(zip(*chg))[1])

	#Cubic-spline interpolation
	cs = CubicSpline(x,y,bc_type='periodic')
	xnew = np.linspace(0, x[-1], grid)
	ynew = list(cs(xnew))

	# calculate ion contribution to the dipole
	M_i, M_e_new = 0, 0
	for e in ion:
		M_i = M_i + e[1]

	# calculate electron contribution to the dipole for interpolated data
	for xx, yy in zip(xnew, ynew):
		M_e_new = M_e_new - yy*xnew[1]
	M_e_new  = M_e_new/ratio

	# total dipole!
	M_tot_new = M_i + M_e_new
	chg.pop(-1)
	return [M_i, M_e_new, M_tot_new]

def getDipole(chg, ion, shift, grid, vol, ratio):
	# add one more data for periodic boundary condition
	chg.append([chg[-1][0] + chg[1][0], chg[0][1]])

	# set x and y
	x = np.asarray(list(zip(*chg))[0])
	y = np.asarray(list(zip(*chg))[1])

	#Cubic-spline interpolation
	cs = CubicSpline(x,y,bc_type='periodic')
	xnew = np.linspace(0, x[-1], grid)
	ynew = list(cs(xnew))
	# print(ynew[0])
	# print(ynew[-2])
	ynew.pop(-1)
	ynew = ynew[-shift:] + ynew[:-shift]
	ynew.append(ynew[0])
	# print(ynew[0])
	# print(ynew[-2])
	# print(len(yyy))
	# charge nomalization regarding added data above
	ynew[0] = ynew[0]/2
	ynew[-1] = ynew[-1]/2

	# calculate ion contribution to the dipole
	D_i, D_e_new = 0, 0
	scaling = 1/(vol*A**3)*q_e*q_e/(2*eps_0)/q_e
	for e in ion:
		D_i = D_i + (e[0] + shift * xnew[1]) * A * e[1]
	D_i = D_i*scaling*A

	# calculate electron contribution to the dipole for interpolated data
	for xx, yy in zip(xnew, ynew):
		D_e_new = D_e_new - xx*A*yy*xnew[1]
	D_e_new = D_e_new*scaling*A/ratio

	# total dipole!
	D_tot_new = D_i + D_e_new
	chg.pop(-1)
	return [D_i, D_e_new, D_tot_new]

def getQuad(chg, ion, shift, grid, vol, ratio):
	# add one more data for periodic boundary condition
	chg.append([chg[-1][0] + chg[1][0], chg[0][1]])

	# set x and y
	x = np.asarray(list(zip(*chg))[0])
	y = np.asarray(list(zip(*chg))[1])

	#Cubic-spline interpolation
	cs = CubicSpline(x,y,bc_type='periodic')
	xnew = np.linspace(0, x[-1], grid)
	ynew = list(cs(xnew))
	ynew.pop(-1)
	ynew = ynew[-shift:] + ynew[:-shift]
	ynew.append(ynew[0])
	# print(len(yyy))
	# charge nomalization regarding added data above
	ynew[0] = ynew[0]/2
	ynew[-1] = ynew[-1]/2

	# calculate ion contribution to the dipole
	Q_i, Q_e_new = 0, 0
	scaling = 1/(vol*A**3)*q_e*q_e/(2*eps_0)/q_e
	for e in ion:
		Q_i = Q_i + ((e[0] + shift * xnew[1]) * A)**2 * e[1]
	Q_i = Q_i*scaling

	# calculate electron contribution to the dipole for interpolated data
	for xx, yy in zip(xnew, ynew):
		Q_e_new = Q_e_new - (xx*A)**2*yy*xnew[1]
	Q_e_new = Q_e_new*scaling/ratio

	# total dipole!
	Q_tot_new = Q_i + Q_e_new
	chg.pop(-1)
	return [Q_i, Q_e_new, Q_tot_new]

def plotChg(chg, shift, grid):
	# add one more data for periodic boundary condition
	chg.append([chg[-1][0] + chg[1][0], chg[0][1]])

	# set x and y
	x = np.asarray(list(zip(*chg))[0])
	y = np.asarray(list(zip(*chg))[1])

	#Cubic-spline interpolation
	cs = CubicSpline(x,y,bc_type='periodic')
	xnew = np.linspace(0, x[-1], grid)
	ynew = list(cs(xnew))
	ynew.pop(-1)
	yshift = ynew[-shift:] + ynew[:-shift]
	yshift.append(yshift[0])
	ynew.append(ynew[0])

	# charge nomalization regarding added data above
	ynew[0] = ynew[0]/2
	ynew[-1] = ynew[-1]/2
	yshift[0] = yshift[0]/2
	yshift[-1] = yshift[-1]/2

	plt.figure()
	# plt.plot(x, y, 'x', xnew, ynew, xnew, np.sin(xnew), x, y, 'b')
	plt.plot(x, y, 'x', xnew, ynew, xnew,yshift)
	plt.legend(['raw', 'Cubic-spline', 'Cubit-spline + shift'])
	# plt.xlim([1.3,1.55])
	# plt.ylim([28,36])
	# plt.axis([-0.05, 0.1, 1.8, 2])
	plt.title('Cubic-spline interpolation')
	plt.show()
	chg.pop(-1)
