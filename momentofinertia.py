import sys
import numpy as np
from matplotlib import pyplot as plt
import math
import astropy as astropy
from astropy import units as u
from astropy.coordinates import SkyCoord

def readin(raw_data):
    raw_list = []
    for line in raw_data:
        ls = line.split(',')
        ls[-1] = ls[-1][:-1]  #removes newline characters at the end
        for i in range(len(ls)):
            try:
                ls[i]  = float(ls[i])
            except:
                pass
        raw_list.append(ls)

    new_list = []

    ### TRANSPOSING

    for j in range(len(raw_list[0])):
        new_list.append([])
        
    for i in range(1, len(raw_list)):
        for j in range(len(raw_list[i])):
            new_list[j].append(raw_list[i][j])

    return new_list

def toXY(name, coords, bcgcoords): 
    """
    Converts RA and dec to x and y by using BCG position as center for tangent approximation
    """
    coordsxy = []

    for i in range(len(coords)):
        ra = coords[i][0]
        dec = coords[i][1]

        comcoords = None
        for j in bcgcoords:
            if name == j[0]:
                comcoords = [j[1], j[2]]
        if comcoords == None:
            comcoords = [0,0]

        ranew = ra - comcoords[0]
        decnew = dec - comcoords[1]

        coordsxy.append([])
        
        x = ranew*math.cos(math.radians(comcoords[1]))
        y = decnew

        coordsxy[i].append(x)
        coordsxy[i].append(y)

    if comcoords == [0,0]:
        print("No BCG for cluster below so likely inaccurate")
    return coordsxy

def parsebcgs(file):
    raw_data = open(file)
    raw_data = raw_data.read()
    lines = raw_data.split("\n")

    data = []
    for i in range(len(lines)):
        elem = lines[i].split()
        data.append(elem)

    if data[-1] == []:
        data = data[:-1]

    bcgcoords = []
    for i in range(len(data)):
        bcgcoords.append([])
        bcgcoords[i].append(data[i][0].lower())
        if bcgcoords[i][0][-1] == '.':
            bcgcoords[i][0] = bcgcoords[i][0][:-1]

        if data[i][3][-1] == '.':
            data[i][3] = data[i][3][:-1]

        if data[i][6][-1] == '.':
            data[i][6] = data[i][6][:-1]


        bcgcoords[i].append(str(data[i][1])+ " " +str(data[i][2]) + " " +str(data[i][3])+ " " + str(data[i][4]) + " " + str(data[i][5]) + " " + str(data[i][6]))

    for i in range(len(bcgcoords)):
        c = SkyCoord(bcgcoords[i][1], unit=(u.hourangle, u.deg))
        bcgcoords[i][1] = c.ra.degree
        bcgcoords[i].append(c.dec.degree)

    return bcgcoords


def momofinertia(coordsxy, weights, name):
    """
    Computes moment of inertia tensor as an nparray, assuming everything is in XY plane
    """    
    center = [0,0]

    ##### ADDING THIS ONLY TO TEST .DAT EXAMPLES#####
        
    #xcenter, ycenter = 6992.611, 2510.27 #A2261

    #xcenter, ycenter = 4599.335, 5442.433 #A381

    xcenter, ycenter = center[0], center[1]         

    #making coordsxy be relative to the center of mass

    for i in range(len(coordsxy)):
        coordsxy[i][0] = coordsxy[i][0]-xcenter
        coordsxy[i][1] = coordsxy[i][1]-ycenter

    momofinertia = np.zeros((3,3), float)
    
    #calculating Ixx:

    Ixx = 0
    for i in range(len(coordsxy)):
        Ixx += weights[i]*coordsxy[i][1]**2
    momofinertia[0][0] = Ixx

    #calculating Iyy: 
    Iyy = 0
    for i in range(len(coordsxy)):
        Iyy += weights[i]*coordsxy[i][0]**2
    momofinertia[1][1] = Iyy


    # calculating Izz:
    Izz = 0
    for i in range(len(coordsxy)):
            Izz += weights[i]*(coordsxy[i][0]**2 + coordsxy[i][1]**2)

    momofinertia[2][2] = Izz

    # calculating Ixy:

    Ixy = 0
    for i in range(len(coordsxy)):
        Ixy += weights[i]*coordsxy[i][0]*coordsxy[i][1]

    momofinertia[0][1] = -Ixy
    momofinertia[1][0] = -Ixy    

    return momofinertia

def getaxis(momofinertia, name, printinfo=True):

    """
    Calculates eigenvalues and eigenvectors of moment of inertia, returns angle
    """
    
    eig = np.linalg.eigh(momofinertia)
    
    eigvals = list(eig[0])
    eigvecs = np.transpose(eig[1])
    mineigval = eigvals.index(min(eigvals))
    principalaxis = eigvecs[mineigval]
    angle = math.degrees(math.atan2(principalaxis[1], principalaxis[0]))

    ###changing angle to match convention:

    newangle = -(90+angle)
    if newangle < -90:
        newangle = 180 + newangle
    if newangle > 90:
        newangle = 180 - newangle

    if printinfo:
        print(name)
        print("moment of inertia tensor: \n" + str(momofinertia))
        print("eigenvalues:" + str(eig[0]))
        print("eigenvectors: \n" + str(eig[1]))
        print("principal axis:" + str(principalaxis))
        print("angle:" + str(newangle))
    
    return eig, angle, newangle

def plot(coordsxy, angle, name):

    plt.subplots()
    x = []
    y = []

    for i in coordsxy:
        x.append(i[0])
        y.append(i[1])

    plt.plot(x, y, 'o', color='black')

    # adding axis segment
    xlim = plt.gca().get_xlim()
    ylim = plt.gca().get_ylim()
    
    plt.gca().set_xlim(xlim)
    plt.gca().set_ylim(ylim)

    xmiddle = (xlim[0] + xlim[1])/2.
    ymiddle = (ylim[0] + ylim[1])/2.
    
    length = 0.05

    endx = length * math.cos(math.radians(angle))
    endy = length * math.sin(math.radians(angle))
    plt.plot((xmiddle, xmiddle+endx), (ymiddle, ymiddle+endy), linestyle="-", color="red", marker='')
    plt.plot((xmiddle, xmiddle-endx), (ymiddle, ymiddle-endy), linestyle="-", color="red", marker='')
    plt.axis("equal")
    plt.gca().invert_xaxis()

    #finding angle according to convention
    newangle = -(90+angle)
    if newangle < -90:
        newangle = 180 + newangle
    if newangle > 90:
        newangle = 180 - newangle
    
    plt.title(name + " angle: " + str(newangle))
    plt.show()

def calculate(coords=None, name=None, plotgals=True, radec=True, printinfo=False):

    if coords==None and name==None: #for file entry mode
        rawdata1 = open(sys.argv[1])
        data = readin(rawdata1)

        #getting name for later
        fullname = str(sys.argv[1])
        allparts = fullname.split("/")
        name = (allparts[-1].split(".")[0]).lower()
        
        #print(name) #for debugging

        #print(data) #for debugging
        
        ra = data[2]
        dec = data[3]
        coords = []

        for i in range(len(ra)):
            coords.append([])
            coords[i].append(ra[i])
            coords[i].append(dec[i])
    else:
        coords = coords
        name = name

    ### Setting all weights equal to 1
    weights = []
    for i in range(len(coords)):
        weights.append(1.)

    ### Getting brightest galaxy coordinates

    bcgcoords = parsebcgs("/home/elenarom/galaxy_alignment/bcgs.txt")

    ### Converting to XY, or not if not necessary
    if radec:
        coordsxy = toXY(name, coords, bcgcoords)
    else:
        coordsxy = coords


    

    momofinertia1 = momofinertia(coordsxy, weights, name)

    angle = getaxis(momofinertia1, name, printinfo)[1]

    newangle = getaxis(momofinertia1, name, printinfo)[2]

    if plotgals:
        plot(coordsxy, angle, name)

    print(name + "," + str(newangle))


calculate(plotgals=False, printinfo=False, radec=True)



##### TEST SET #####

# testcoordsxy1 = [[0,0], [1,1], [2,2], [3,3], [4,4]]
# testcoordsxy = [[0,0], [2,4], [2,0], [0,4]]

# momofinertiatest = momofinertia(testcoordsxy, False)
# angle = getaxis(momofinertiatest, "test", True, False)[1]
# plot(testcoordsxy, angle, "test", False)



# raw_data = open("/home/elenarom/galaxy_alignment/a383par.dat")
# raw_data = raw_data.read()
# lines = raw_data.split("\n ")
# data = []
# for i in range(len(lines)):
#     elem = lines[i].split()
#     data.append(elem)

# testcoordsxy = []

# for i in range(len(data)):
#     testcoordsxy.append([])
#     testcoordsxy[i].append(float(data[i][0]))
#     testcoordsxy[i].append(float(data[i][1]))


#calculate(coords=testcoordsxy, name="a383 test", printinfo=True, radec=False)
    


