#! /usr/bin/env python
# coding: utf-8


import numpy as np
#import pandas as pd
import re
#from itertools import cycle
#import  operator
from collections import OrderedDict
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pylab

import cPickle as pickle #needed while saving data
import time
import whole_brain_neuron_name as nn
import os
import random
start = time.time()
xyz = {}
class Scalar:
    value = 0.5
    sigma = 3.0


def radius(x, y, z, xm, ym, zm): #distance between two points
    r = np.sqrt((x-xm)**2+(y-ym)**2+(z-zm)
                **2)
    return r

def prob(r, sigma): #probability formula
    p = np.exp(-(r/sigma)**2/2)/(np.sqrt(2*(np.pi)**3)*(sigma)**3)
    return p

def toScalar(x): #convert to resolution
    if round((x%Scalar.value)/Scalar.value) == 1:
        x = x-x%Scalar.value+Scalar.value
    else:
        x = x-x%Scalar.value
    return float(x)  # actually float(x) already present the same value as float("%.3f" %x)


def name_to_index(name): #convet name to int
    return nn.namelist.index(name)

def pos_to_int(x): #convert to int(x*10)
    return int(round(10*x))

def theta_to_int(x): #convert to int(x*100)
    return int(round(100*x))


countpoints = 0
allpoints = OrderedDict()
exist = set()  # ***check if exist is needed to be outside the for loop
xyz=[]
time = 0
swc = " "
count=0

for dirPath, dirNames, fileNames in os.walk("bundle_space_100_new"):
   for f in fileNames:
        print os.path.join(dirPath, f)

url=fileNames


x2=[]
y2=[]
z2=[]
check=0

regex_2 = re.compile(
    r"([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\d]+)[\s]+([\d]+)")
index=0
previous_id=-1
previous_branch=-1
pathway='bundle_space_100_fig_100/'
fig = plt.figure()
ax = fig.gca(projection='3d')
# spxn1yp3zp2.txt
for space in url:
    if index==3:
        break
    index=index+1
    print space
    position_x=4
    position_y=space.index("y")+2
    position_z=space.index("z")+2
    if space[3]=='n':
        sign_x=-1
    else:
        sign_x=1
    if space[position_y-1]=='n':
        sign_y=-1
    else:
        sign_y=1
    if space[position_z-1]=='n':
        sign_z=-1
    else:
        sign_z=1
    x_range=int(space[position_x:(position_y-2)])
    y_range=int(space[position_y:(position_z-2)])
    z_range=int(space[position_z:-4])
    print "x_range", sign_x, x_range
    print "y_range", sign_y, y_range
    print "z_range", sign_z, z_range
    # ax.set_xlim(sign_x * x_range * 100, sign_x * (x_range + 1) * 100)
    # ax.set_ylim(sign_y * y_range * 100, sign_y * (y_range + 1) * 100)
    # ax.set_zlim(sign_z * z_range * 100, sign_z * (z_range + 1) * 100)
    # ax.set_xlim(0,1000)
    # ax.set_ylim(0,1000)
    # ax.set_zlim(0,1000)

    with open("bundle_space_100_new/"+space) as space_file:
        for line in space_file:
            # print line
            for iter2 in regex_2.finditer(line):
                id=float(iter2.group(4))
                branch_id=float(iter2.group(5))
                if previous_id!=id or branch_id!= previous_branch:
                    r = random.random()
                    g = random.random()
                    b = random.random()
                    plt.plot(x2, y2, z2,color=(r, g, b, 0.5))
                    # print x2
                    # print y2
                    # print z2
                    # print previous_id
                    # print previous_branch
                    previous_id = id
                    previous_branch = branch_id

                    x2=[]
                    y2=[]
                    z2=[]
                    check=1
                x2.append(float(iter2.group(1)))
                y2.append(float(iter2.group(2)))
                z2.append(float(iter2.group(3)))
        r = random.random()
        g = random.random()
        b = random.random()
        # print x2
        # print y2
        # print z2
        # print id
        # print branch_id
        plt.plot(x2, y2, z2, color=(r, g, b, 0.5))
        previous_id = id
        x2 = []
        y2 = []
        z2 = []
        print space[:-4]
        # ax.xlim(sign_x*x_range*100, sign_x*(x_range+1)*100)
        # ax.ylim(sign_y*y_range*100, sign_y*(y_range+1)*100)
        # ax.zlim(sign_z*z_range*100, sign_z*(z_range+1)*100)
        # ax.axis([sign_x*x_range*100, sign_x*(x_range+1)*100, sign_y*y_range*100, sign_y*(y_range+1)*100,sign_z*z_range*100, sign_z*(z_range+1)*100])
        # plt.savefig(pathway+space[:-4]+'.tiff', dpi = 100)
        check=0
        x2=[]
        y2=[]
        z2=[]
plt.show()
plt.close()
fig = plt.figure()
ax = fig.gca(projection='3d')


