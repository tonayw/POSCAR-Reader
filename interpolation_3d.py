import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import RegularGridInterpolator as RGI
import sys, getopt
import math
from decimal import Decimal
from itertools import islice

def setPOSCARinfo(POSCAR, ionCHG):
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

    # set height
    cosPHI = vol/(np.linalg.norm(np.cross(pVec1, pVec2))*a3)
    h = a3*cosPHI
    print(h)

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
    #    crd[ii] = list(map(float,lines[ii+8].split()))[2]*h
        ii = ii + 1

    ion = [["",""] for e in range(N)]
    if(coord == "Direct"):
        for e, i, cr in zip(ion, ionchg, crd):
            e[0] = cr
            e[1] = i
    else:
        print("Cartesian coordinate not supported!!")
    return h, vol, ion

# when CHGCAR file is large, use this method
def setCHGLarge(CHGCAR,ionCHG):
    with open(CHGCAR) as f:
        # save title 
        for line in islice(f, 1):
            title = line
        # print("header: " % (title))
        # set scaling factor
        for line in islice(f, 1):
            sf = float(list(map(float,line.split()))[0])
        i = 0

        # set primitive vectors
        pVec = [[None for _ in range(3)] for _ in range(3)]
        for line in islice(f, 3):
            pVec[i] = list(map(float,line.split()))
            i = i + 1
        # multiply scaling factor to the lattice vectors
        pVec[0] = [e * sf for e in pVec[0]]
        pVec[1] = [e * sf for e in pVec[1]]
        pVec[2] = [e * sf for e in pVec[2]]

        # set lattice parameters
        a = [0 for i in range(3)]
        a[0] = np.linalg.norm(pVec[0])
        a[1] = np.linalg.norm(pVec[1])
        a[2] = np.linalg.norm(pVec[2])
        print("lattice parameters: (%.4f, %.4f, %.4f)" % (a[0], a[1], a[2]))
        # set atoms
        tmp = []
        for line in islice(f, 2):
            tmp.append(line)
        atoms = tmp[0].split()      # array of atom symbols
        atomN = list(map(int,tmp[1].split()))    # array of # of atoms per type
        print("atom info.:",atoms, atomN)
        # print(atomN)
        # set # of atoms: N
        N = 0
        for e in atomN:
            N = N + e
        # print(N)

        # set atom ion charges for each component
        ionchg = [[""] for e in range(N)]
        ii, jj = 0, 0
        for atom in atomN:
            for e in range(atom):
                ionchg[ii] = ionCHG[jj]
                ii = ii + 1
            jj = jj + 1

        # set cell volume
        vol = np.dot(pVec[2], np.cross(pVec[0], pVec[1]))
        print("volume: ", vol)
        # set height
        cosPHI = vol/(np.linalg.norm(np.cross(pVec[0], pVec[1]))*a[2])
        h = a[2]*cosPHI
        print("h:", h)
        # set Direct or Cart
        for line in islice(f, 1):
            if(line[0] == "D" or line[0] == "d"):
                coord = "Direct"
            elif(line[0] == "D" or line[0] == "d"):
                coord = "Cart"
            else:
                print("POSCAR is neither Direct nor Cart. SOMETHING IS WRONG!!")
                coord = "SOMETHING WRONG"

        # set Atomic positions
        atomCrd = []
        for line in islice(f,N):
            line = list(map(float,line.split()))
            atomCrd.append(line)

        # set Nz = grid[0], Ny = grid[1], Nz = grid[2]
        for line in islice(f,2):
            grid = list(map(int,line.split()))

        # total # of grids
        gridN = grid[0]*grid[1]*grid[2]
        print("grids: ", grid)

        # write CHGCAR
        # out = open('CHGCAR_test', 'w')
        # out.write("%s" % (title))
        # out.write("   %.14f\n" % (1))
        # for e in range(3):
        #     out.write("     %.6f    %.6f    %.6f\n" % (pVec[e][0],pVec[e][1],pVec[e][2]))
        # for at in atoms:
        #     out.write("   %s" % at)
        # out.write("\n")
        # for at in atomN:
        #     out.write("   %d" % at)
        # out.write("\n")
        # out.write("%s\n" % (coord))

        # for at in atomCrd:
        #     out.write("  %.6f  %.6f  %.6f\n" % (at[0],at[1],at[2]))
        # out.write("\n   %d   %d   %d\n" % (grid[0],grid[1],grid[2]))

        # set CHG data
        chgsum = 0
        chgdat = [[[0 for k in range(grid[2])] for j in range(grid[1])] for i in range(grid[0])]
        i = 0
        for line in islice(f,math.ceil(gridN/5.0)):
            tmp = list(map(float,line.split())) 
            j = 0
            for t in tmp:
                nn = i*5 + j
                g3 = math.floor(nn/(grid[0]*grid[1]))
                gg = nn%(grid[0]*grid[1])
                g2 = math.floor(gg/grid[0])
                g1 = nn%grid[0]
                chgdat[g1][g2][g3] = t
                j = j + 1
            i = i + 1
            # for t in tmp: out.write(" %0.010E" % (t))
        print("total charge: %.5f" % (chgsum/grid[0]/grid[1]/grid[2]))
        return [chgdat, grid, pVec, a, vol]

def writCHG():
    f = open('CHGCAR_test', 'w')
    
    f.write("   %d   %d   %d \n" % (grid[0],grid[1],grid[2]))
    f.close()

# add additional CHG layer at boundaries for proper interpolation
def addCHGlayer(chgdat,grid):
    chg = [[[0 for k in range(grid[2]+1)] for j in range(grid[1]+1)] for i in range(grid[0]+1)]
    for i in range(len(chg[0][0])):
        for j in range(len(chg[0])):
            for k in range(len(chg)):
                if i < grid[0] and j < grid[1] and k < grid[2]:
                    chg[i][j][k] = chgdat[i][j][k]
                else:
                    if i == grid[0] and j <  grid[1] and k <  grid[2]:
                        chg[i][j][k] = chgdat[0][j][k]
                    if i <  grid[0] and j == grid[1] and k <  grid[2]:
                        chg[i][j][k] = chgdat[i][0][k]
                    if i <  grid[0] and j <  grid[1] and k == grid[2]:
                        chg[i][j][k] = chgdat[i][j][0]

                    if i <  grid[0] and j == grid[1] and k == grid[2]:
                        chg[i][j][k] = chgdat[i][0][0]
                    if i == grid[0] and j <  grid[1] and k == grid[2]:
                        chg[i][j][k] = chgdat[0][j][0]
                    if i == grid[0] and j == grid[1] and k <  grid[2]:
                        chg[i][j][k] = chgdat[0][0][k]

                    if i == grid[0] and j == grid[1] and k == grid[2]:
                        chg[i][j][k] = chgdat[0][0][0]  
    return chg

# def shellChgDen(chgdat,r):
def shellChgDen(r,chgdat,grid,pVec,a):
    chg = addCHGlayer(chgdat,grid)
    x = np.linspace(0,1,grid[0]+1)
    y = np.linspace(0,1,grid[1]+1)
    z = np.linspace(0,1,grid[2]+1)
    # print(x)
    func = RGI((x,y,z),chg)
    # pts = (1/grid[0]*a[0],0,0)
    # pts = (0,0,0)
    # print(func(pts))
    thGrid = int(r*400)
    if thGrid < 200: thGrid = 200
    shellChgSum = 0
    pointN = 0
    for i in range(thGrid):
        th = math.pi*i/thGrid
        phGrid = int(round(thGrid*2*math.sin(th)))
        for j in range(phGrid):
            ph = 2*math.pi*j/(phGrid)

            x = r * math.sin(th) * math.cos(ph)
            y = r * math.sin(th) * math.sin(ph)
            z = r * math.cos(th)

            vec = [x,y,z]

            # set r center vector
            r1 = [e * 0.5 for e in pVec[0]]
            r2 = [e * 0.5 for e in pVec[1]]
            r3 = [e * 0.5 for e in pVec[2]]
            rcenter = list(map(sum, zip(r1,r2,r3)))

            v_c = list(map(sum, zip(vec,rcenter)))

            M = np.array([[pVec[0][0], pVec[1][0], pVec[2][0]],
                            [pVec[0][1], pVec[1][1], pVec[2][1]],
                            [pVec[0][2], pVec[1][2], pVec[2][2]]])
            inverseM = np.linalg.inv(M)

            r_c = inverseM.dot(v_c)

            # print(M)
            # print(r_c)
            # [i - j for i, j in zip(a, b)]
            # vecFromCenter = list(map(sum, zip(rcenter,vec)))
            # print(vecFromCenter)
            r_c[0] = r_c[0]%1.0
            r_c[1] = r_c[1]%1.0
            r_c[2] = r_c[2]%1.0

            pts = (r_c[0],r_c[1],r_c[2])
            # print(func(pts))
            shellChgSum += func(pts)
            pointN += 1
            # print(ph)
    rho = 0
    print(r, shellChgSum/pointN)
    # print(shellChgSum/pointN)

    return rho
# data = setPOSCARinfo("POSCAR",[4])
# print(data[0], data[1], data[2])
# chg = setSphCHG("CHGCAR_sum",[6],1)

def getPointCHG(point,chgdat,grid,pVec,a):
    chg = addCHGlayer(chgdat,grid)
    x = np.linspace(0,1,grid[0]+1)
    y = np.linspace(0,1,grid[1]+1)
    z = np.linspace(0,1,grid[2]+1)
    # print(x)
    func = RGI((x,y,z),chg)
    # pts = (1/grid[0]*a[0],0,0)
    # pts = (0,0,0)
    # print(func(pts))

    pointCHG = func(point)

    return pointCHG
# multiple points
def getPointsCHG(lis_points,chgdat,grid,pVec,a):
    chg = addCHGlayer(chgdat,grid)
    x = np.linspace(0,1,grid[0]+1)
    y = np.linspace(0,1,grid[1]+1)
    z = np.linspace(0,1,grid[2]+1)
    # print(x)
    func = RGI((x,y,z),chg)
    # pts = (1/grid[0]*a[0],0,0)
    # pts = (0,0,0)
    # print(func(pts))
    lis = []
    for p in lis_points:
        lis.append(pointCHG = func(p))
    return lis

if __name__ == "__main__":
    material = 'c-C'
    folder = 'data'
    #fileName = './'+folder+'/'+material+'_CHGCARdiff'    # planar averaged charge density
    fileName = './'+'c-BN_CHGCARdiff_light'
    ionCHG = [3,5] # valence electrons   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!               WHAT IS THIS???????????

    data = setCHGLarge(fileName,ionCHG)
    #data = setPOSCARinfo(fileName,ionCHG)

    # set a point in a direct coordinate (0 to 1)
    x1_d, x2_d, x3_d = 0, 0, 0

    point = (x1_d%1.0, x2_d%1.0, x3_d%1.0) 
    pointCHG = getPointCHG(point,data[0],data[1],data[2],data[3])
    print(type(pointCHG))
    print("a given point: ", point)
    print("CHG at a given point: %.5f" % (pointCHG/data[4]))
