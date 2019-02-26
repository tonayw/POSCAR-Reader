import numpy as np

##The filereader returns a list of information about the unit cell in the following form
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
    infolist.append(f.read().split())
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

def average_charge_along_axis(axis,coord,filename):
    total = 0
    two_vectors = []
    #this is a list of arrays
    for i in range(3):
        if i==axis:
            continue
        else:
            two_vectors.append(coord[i])
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
    print(infolist)

    #a_charge_i = average_charge_along_axis(axis,coord,'initial_data.txt')
    #a_charge_f = average_charge_along_axis(axis,coord,'final_data.txt')