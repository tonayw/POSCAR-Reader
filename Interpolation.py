'''
Interpolation:
Import from 3Dinterpolation

Interpolation along direction between two atoms

Atom1 position, Atom 2 position ---> gives the vector pointing from atom 1 to atom 2 -----> find the unit vector


N --> accuracy of data point along this vector

point coordinate (along that specific vector) = vector * n/N

using vector addition

(need more details ----------> discussion on meeting)

'''

from interpolation_3d import *




def transfer_to_vector(tail,head):
    vector = []
    vector.append(tail[0]-head[0])
    vector.append(tail[1]-head[1])
    vector.append(tail[2]-head[2])
    return vector

'''The filereader returns a list of information about the unit cell in the following form:
      [scaling factor, first axis, second axis, third axis, [first atom position, second atom position...],
      [grid dimensions],[data points]]
'''
def CHGCAR_filereader(filename):
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
    atomnumlist = f.readline().split()
    print(atomnumlist)
    for i in range(len(atomnumlist)):
        print(type(atomnumlist[i]))
        atomnumlist[i] = int(atomnumlist[i])
    print(atomnumlist)
    atomnumber = sum(atomnumlist)
    print(atomnumber)
    next(f)
    atomlist = []
    line = f.readline()
    for i in range(atomnumber):
        smallatomlist = line.split()
        for j in range(3):
            smallatomlist[j] = float(smallatomlist[j])
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
    
def distance(point1,point2):
    D = np.sqrt((point1[0]-point2[0])**2+(point1[1]-point2[1])**2+(point1[2]-point2[2])**2)
    return D

def interpoate(point1,point2,n,data):
    p1p2 = transfer_to_vector(point2,point1) #p1p2 = [x2-x1,y2-y1,z2-z1]
    #print(point1,point2)
    d = distance(point1,point2)
    chglist = []
    for i in range(n):
        coord0 = point1[0]+i/n*p1p2[0]
        coord1 = point1[1]+i/n*p1p2[1]
        coord2 = point1[2]+i/n*p1p2[2]
        if i == n-1:
            print(coord0, coord1, coord2)
        
        point = (coord0%1.0,coord1%1.0,coord2%1.0)
        pointcharge = getPointCHG(point,data[0],data[1],data[2],data[3])
        chglist.append([point,pointcharge/data[4]])
    return chglist
        
        
filename = "CHGCAR_sum"
ionchg = [3,5]
x = CHGCAR_filereader(filename)
data = setCHGLarge(filename,ionchg)
print(data[4])

accuracy = 500
atomlist = x[4]
print(x[4])

interpolation_result = []
n = 0


for i in range(len(atomlist)):
    for j in range(n,len(atomlist)):
        if (i!=j):
            interpolation_result.append(interpoate(atomlist[i],atomlist[j],accuracy,data))
    n+=1    
print(interpolation_result)



for i in range(len(interpolation_result)):
    plotting_array = np.zeros(len(interpolation_result[i]))
    print(i)
    for j in range(len(interpolation_result[i])):
        print("point:",j,"           charge: ",interpolation_result[i][j][1])
'''
    plt.figure()
    plt.plot(newlist[0],newlist[1])
    plt.title("Charge Density")
    plt.xlabel("r")
    plt.ylabel("E")
    plt.show()
'''