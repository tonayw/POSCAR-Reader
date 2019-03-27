import numpy as np
import matplotlib
import sys
import matplotlib.pyplot as plt

'''This returns a list that contains all information from the POSCAR file
The list is of the following form:
[scaling factor,[first axis,second axis, third axis](multiplied by scaling factor),
[[first element, number of atoms],[second element, number of atoms]...],
[coordinate of first atom, seocnd atom...](listed in sequence of the atoms above)]

'''

def diagnal(infolist):
    diagnalv = [0.0,0.0,0.0]
    for i in range(3):
        for j in range(3):
            diagnalv[i]+= infolist[1][i][j]
    return np.sqrt(diagnalv[0]**2+diagnalv[1]**2+diagnalv[2]**2)
    
def filereader(filename):
    infolist = []
    f = open(filename)
    line =f.readline()
    print("\n")
    print("name:"+line)
    print("scaling factor:")
    line = f.readline()
    print(line)
    scalingfactor = float(line)
    infolist.append(scalingfactor)
    coordinatelist = []
    for i in range (3):
        vectorlist = f.readline().split()#save the vector to a list
        for j in range(3):
            vectorlist[j] = float(vectorlist[j])*scalingfactor
        coordinatelist.append(vectorlist)
    infolist.append(coordinatelist)
    atoms = []#we use this list to save the information of the atoms
    names = f.readline().split()  #the list of atom names
    num_atoms = f.readline().split()#the list of the associated atom numbers
    for i in range(len(names)):
        singleatom=[]
        singleatom.append(names[i])
        singleatom.append(int(num_atoms[i]))
        atoms.append(singleatom)
    infolist.append(atoms)
    next(f)
    positionlist = []
    line = f.readline()
    while line:
        linelist = line.split()#linelist contains three original coordinates
        #from the file that gives the position of the atom
        for i in range(len(linelist)):
            linelist[i] = float(linelist[i])
        atom_coord = []#this stores the cartesian coordinates of single atom
        for i in range(3):
            final_coord = 0 #reset the accumulating coordinate every time           
            for j in range(3):
                final_coord = final_coord+linelist[i]*coordinatelist[j][i]
            atom_coord.append(final_coord)
        positionlist.append(atom_coord) 
        line = f.readline()
    infolist.append(positionlist)
    return infolist

def coordplus(coord1,coord2):
    return [coord1[0]+coord2[0],coord1[1]+coord2[1],coord1[2]+coord2[2]]

def coordminus(coord1,coord2):
    return [coord1[0]-coord2[0],coord1[1]-coord2[1],coord1[2]-coord2[2]]

def distance(point1,point2):
    D = np.sqrt((point1[0]-point2[0])**2+(point1[1]-point2[1])**2+(point1[2]-point2[2])**2)
    return D

def bond_calculator(infolist,lower_bound,upper_bound):
    bondlist = []
    total_bond = 0
    criteria = upper_bound-lower_bound
    #for i in range(len(infolist[3])):
        #bondlist.append([])
    #counting total bond of an atom
    for i in range(len(infolist[3])):
        bond_tot = 0
        for j in range(len(infolist[3])):
            if i == j:
                continue
            else:
                otheratom = infolist[3][j].copy()
                adjusted_otheratom = [0,0,0]
                for x in range(-1,2):
                    for y in range(-1,2):
                        for z in range(-1,2): 
                            adjusted_otheratom[0] = otheratom[0]+infolist[1][0][0]*x+infolist[1][1][0]*y+infolist[1][2][0]*z
                            adjusted_otheratom[1] = otheratom[1]+infolist[1][0][1]*x+infolist[1][1][1]*y+infolist[1][2][1]*z
                            adjusted_otheratom[2] = otheratom[2]+infolist[1][0][2]*x+infolist[1][1][2]*y+infolist[1][2][2]*z
                            if lower_bound<=distance(infolist[3][i],adjusted_otheratom) and upper_bound>=distance(infolist[3][i],adjusted_otheratom):
                                bond_tot+=1
                            else:
                                continue
                        
                    
        #efficient way
        total_bond+=bond_tot
        bondlist.append(bond_tot)
    bondlist.append(total_bond/2)
    return bondlist

def histogram_plotter(infolist,bondlist):
    bondlist_iterator=0 
    histogram_data = []
    for i in range(len(infolist[2])):
        num_bins = 7
        single_atom_bond = []
        for j in range(int(infolist[2][i][1])):
            single_atom_bond.append(bondlist[bondlist_iterator])
            bondlist_iterator+=1
        histogram_data.append(single_atom_bond)
        n, bins, patches = plt.hist(single_atom_bond, num_bins, facecolor='blue', alpha=0.5)
        plt.show()
    return 0

if __name__ == '__main__':
    x = filereader("POSCAR-GST")
    #x = filereader("POSCAR-BaTiO3")
    #x = filereader("POSCAR-Si")
    print(x)
    #y = bond_calculator(x,2,2.05)
    #y = bond_calculator(x,2.3,2.4)#Si
    y = bond_calculator(x,2.5,3.6)
    print("Total number of bonds is: ")
    print(y[-1])
    histogram_plotter(x,y)
    