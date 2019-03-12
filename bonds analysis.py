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
def filereader(filename):
    infolist = []
    f = open(filename)
    line =f.readline()
    print("\n")
    print("name:"+line)
    print("scaling factor:")
    line = f.readline()
    scalingfactor = float(line)
    infolist.append(scalingfactor)
    print(line)
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
        singleatom.append(num_atoms[i])
        atoms.append(singleatom)
    infolist.append(atoms)
    next(f)
    positionlist = []
    line = f.readline()
    print(coordinatelist)
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

def distance(point1,point2):
    D = np.sqrt((point1[0]-point2[0])**2+(point1[1]-point2[1])**2+(point1[2]-point2[2])**2)
    return D

#def bond_calculator(info_list,lower_bound,upper_bound):
if __name__ == '__main__':
    x = filereader("POSCAR-BaTiO3")
    print(x)
    