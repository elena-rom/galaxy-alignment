import numpy as np
from matplotlib import pyplot as plt
import math
import astropy as astropy
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy import stats as stats

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

    return raw_list

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

    for i in range(len(data)):
        for j in range(7,11):
            bcgcoords[i].append(data[i][j])
            if bcgcoords[i][j-4][-1] == '.':
                bcgcoords[i][j-4] = bcgcoords[i][j-4][:-1]
            bcgcoords[i][j-4] = float(bcgcoords[i][j-4])

    return bcgcoords

inertiaclusterangles = readin(open("result.csv"))

"""
inertiaclusterangles data:
0: name
1: cluster moment of inertia angle
"""

bcgangles = readin(open("bcg_angles.csv"))[1:]

"""
bcgangles data:
0: name
1: bcg angle
2: bcg angle error
3: bcg ellipticity
"""

### making all names lowercase for parsing purposes

for i in bcgangles:
    i[0] = i[0].lower()

bcgdata = parsebcgs("bcgs.txt")
"""
bcgdata data:
0: name
1: bcg ra
2: bcg dec
3: cluster mean LOS (line-of-sight) velocity
4: cluster velocity dispersion
5: bcg LOS velocity
6: significance of difference bw cluster and bcg velocity (see email)
"""


### note: each of these has a slightly different collection of galaxy clusters, so it's necessary to collect each by the name string


anglecomparison = [] ### list of [cluster name, momofinertia angle, bcg angle, diff significance]

for i in range(len(bcgangles)):
    name = bcgangles[i][0]
    for j in range(len(bcgdata)):
        for k in range(len(inertiaclusterangles)):
            if name == bcgdata[j][0] and name == inertiaclusterangles[k][0]: #checking that name matches for all 3 sets
                anglecomparison.append([name, inertiaclusterangles[k][1], bcgangles[i][1], bcgdata[j][6]])


### fixing sign convention for bcg angle
# for i in anglecomparison:
#     angle = i[2]
#     newangle = -(90+angle)
#     if newangle < -90:
#         newangle = 180 + newangle
#     if newangle > 90:
#         newangle = 180 - newangle
#     i[2] = newangle


### setting up to make a plot

angle_differences = []
velocity_diff_significances = []

for i in anglecomparison:
    angle = abs(i[2]-i[1])
    if angle > 90:
        if i[2] < 0:
            newcomparison = 180 + i[2]
            angle = abs(newcomparison - i[1])
        if i[1] < 0: 
            newcomparison = 180 + i[1]
            angle = abs(newcomparison - i[2])

    angle_differences.append(angle)
    velocity_diff_significances.append(abs(i[3]))

plt.plot(velocity_diff_significances, angle_differences, 'o', color='black')
plt.gca().set_xlim(-0.2,3)
plt.xticks(np.arange(-0.2, 3, 0.2))
plt.grid(True)
plt.xlabel("BCG/Cluster Velocity Difference")
plt.ylabel("BCG/Cluster Angle Difference")
plt.title("Relationship between Angle and Velocity Differences between BCG and Cluster")
plt.show()

print("N=", len(anglecomparison))

def ks_in_range(data1, data2, lowerbound, upperbound):
    data_in_range = []
    for i in range(len(data1)):
        if data2[i] > lowerbound and data2[i] < upperbound:
            data_in_range.append(data1[i])
    test_stat_range = stats.kstest(data_in_range, 'uniform', args=(0,90))

    print("test stat for range", lowerbound, upperbound, test_stat_range)



test_stat_overall = stats.kstest(angle_differences, 'uniform', args=(0,90))
print("test_stat:", test_stat_overall)

for i in range(0,3):
    ks_in_range(angle_differences, velocity_diff_significances, i*0.2, (i+1)*0.2)












