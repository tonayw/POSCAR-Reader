import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import CubicSpline
import multipoleFunc as mF
import sys, getopt

# Global variables
eps_0 = 8.8541878176E-12
q_e = 1.60217662E-19
A = 1E-10

# option variable
newGrid = 2001

def main(argv):
    # material = 'polar'
    material = 'c-C'
    lattice = 'ZB'        # ZB: Zinc Blende, W: Wurzite, F: FCC, and o: others
    folder = 'data'
    #ionCHG = [10,12,6] 
    ionCHG = [4] # valence electrons

    #fileName = './'+folder+'/'+material+'_chgave.txt'    # planar averaged charge density
    fileName = './'+'c-BN_CHGCARdiff_light_test_version'
    print(fileName)
    POSCAR = './'+folder+'/'+material+'_POSCAR'        

    # latInfo = [lattice, lat_a, lat_c, h, vol, ratio, ion]
    print("planar ave. CHG filename: %s" % fileName)
    print("POSCAR filename:   %s" % POSCAR)

    latInfo = [lattice] + list(mF.setPOSCARinfo(POSCAR, lattice, ionCHG))
    newCHG = getNewCHG(fileName, latInfo)

    # print output
    for e in newCHG:
        print(e[0], e[1])

    # write output
    out = open('./'+folder+'/'+material+'_chgave_interpol.txt', 'w')
    for e in newCHG:
        out.write(str(e[0]) + "\t" + str(e[1]) + "\n")
    out.close()


def getNewCHG(CoreFileName, latInfo):
    cfilename = CoreFileName
    # vfilename = ValFileName    
    lattice = latInfo[0]
    lat_a, lat_c, h, vol = latInfo[1], latInfo[2], latInfo[3], latInfo[4]
    ratio = latInfo[5]
    ion = latInfo[6]
    # print(ion)
    Area = vol/h*A*A

    # print(ion)
    #read file and parse!
    cfile = open(cfilename,"r")
    clines = cfile.read().split("\n")

    #save the data to an array
    cchgdat = [["",""] for e in range(len(clines))]
    
    ii = 0
    for l in clines:
        values = l.split("\t")
        cchgdat[ii][0] = ii*h/len(clines)
        cchgdat[ii][1] = float(values[1])
        ii=ii+1

    chgdat = [cchgdat[i] for i in range(len(cchgdat))]
    chgdat = [[x[0], x[1]] for x in chgdat]

    # calculate multipoles from input data
    tmp = mF.getOriginalMultipoles(chgdat,ion,vol,ratio)
    M_raw, D_raw, Q_raw = tmp[0], tmp[1], tmp[2]

    print("lattice parameter (conventional %s): \n  a = %.5f, c = %.5f" % (lattice, lat_a, lat_c))
    print("cell info: \n  hight = %.5f A, volume: %.5f A^3" % (h,vol))
    print("--------------------------- INPUT CHG INFO -------------------------")
    print(" Charge (ion)   Charge (el)    Dipole (eV/A)   Quadrupole (eV)")
    print("  %+.5f      %+.5f      %+.5f        %.5f " % (M_raw[0],M_raw[1],D_raw[2],Q_raw[2]))
    print("--------------------------------------------------------------------")
    print("Interpolation start!")
    print("  grid: %d -> %d " % (len(chgdat), newGrid))
    # calculate interpolated chgdat
    xnew, ynew = getInterpolation(chgdat,ion,newGrid,ratio)
    newCHG = list(zip(xnew,ynew))
    return newCHG

def getInterpolation(chg,ion,grid,ratio):
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
    return xnew, ynew

if __name__ == "__main__":
    main(sys.argv[1:])