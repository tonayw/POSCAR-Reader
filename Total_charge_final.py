import numpy as np
import matplotlib
import sys
import matplotlib.pyplot as plt

'''The filereader returns a list of information about the unit cell in the following form:
      [scaling factor, first axis, second axis, third axis, [first atom position, second atom position...],
      [grid dimensions],[data points]]
'''
def filereader(filename):
    #infolist contains three vectors and the scaling factor
    infolist = []
    f = open(filename)
    line = f.readline()
    print("name: "+line)
    print("scaling factor:")
    line = f.readline()
    scalingfactor = float(line)
    infolist.append(scalingfactor)
    print(line)
    for i in range (3):
        vectorlist = f.readline().split()#save the vector to a list
        for j in range(3):
            vectorlist[j] = float(vectorlist[j])*scalingfactor
        infolist.append(vectorlist)
    next(f)
    print("atom number:")
    atomnumber = int(f.readline())
    print(atomnumber)
    next(f)
    atomlist = []
    line = f.readline()
    for i in range(atomnumber):
        smallatomlist = line.split()
        for j in range(3):
            smallatomlist[j] = float(smallatomlist[j])*scalingfactor
        atomlist.append(smallatomlist)
        line = f.readline()
    infolist.append(atomlist)
    
    grid = f.readline().split()
    for i in range(3):
        grid[i] = int(grid[i])
    infolist.append(grid)
    print("*************")
    data_list = f.read().split()
    for i in range(len(data_list)):
        data_list[i]=float(data_list[i])
    infolist.append(data_list)
    return infolist
    
def total_charge_calculator(filename):
    
    filename = open(filename)
    total_charge = 0
    for i in range(13):
        next(filename)
    line = filename.readline()
    while line:
        linelist = map(float, line.split())
        for datapt in linelist:
            total_charge+=datapt
        line = filename.readline()
    filename.close()
    return total_charge

def average_charge_along_axis(axis,infolist):
    average_charge_list = []
    tot_chrg = 0
    chg_density = 0
    gridlist = infolist[5].copy()
    axis_grid = gridlist.pop(axis-1)
    print(axis_grid)
    grids_per_plane = gridlist[0]*gridlist[1]
    print("axis length:",np.sqrt(infolist[axis][0]**2+infolist[axis][1]**2+infolist[axis][2]**2))
    unit_length = np.sqrt(infolist[axis][0]**2+infolist[axis][1]**2+infolist[axis][2]**2)/axis_grid
    print(unit_length)
    for i in range(axis_grid):     #three dimensional list
        for j in range(grids_per_plane):
            tot_chrg+=infolist[-1][i*grids_per_plane+j]
        chg_density = tot_chrg/grids_per_plane/(unit_length*axis_grid)
        average_charge_list.append([i*unit_length,chg_density])
        tot_chrg = 0
        chg_density = 0        
    return average_charge_list
            
    
    #this is a list of arrays
    #cross_product = np.cross(two_vectors[0],two_vectors[1])
    #cross_area = np.linalg.norm(cross_product)


if __name__=='__main__':
    #coord is a list of all arrays
    a1 = np.array([0,0.5,0.5])
    a2 = np.array([0.5,0,0.5])
    a3 = np.array([0.5,0.5,0])
    coord = [a1,a2,a3]
    print('\n')
    print("initial total charge:")
    print(total_charge_calculator('CHGCAR-initial2')/32/32/32)
    print('\n')
    print("final total charge:")
    print(total_charge_calculator('final_data.txt')/32/32/32)
    infolist = filereader("CHGCAR-initial - Copy.txt")
    average_charge = average_charge_along_axis(3,infolist)
    newlist =[[],[]]
    for i in range(len(average_charge)):
        print("{:>12}  {:12}".format(average_charge[i][0],average_charge[i][1]))
        newlist[0].append(average_charge[i][0])
        newlist[1].append(average_charge[i][1])
    
    plt.figure()
    plt.plot(newlist[0],newlist[1])
    plt.title("Electric Field")
    plt.xlabel("r")
    plt.ylabel("E")
    plt.show()
    #a_charge_i = average_charge_along_axis(axis,coord,'initial_data.txt')
    #a_charge_f = average_charge_along_axis(axis,coord,'final_data.txt')