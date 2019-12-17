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
    r = np.sqrt((x-xm)**2+(y-ym)**2+(z-zm)**2)
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
regex = re.compile(
    r"([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\d]+)"
)

# ############somata
fig = plt.figure()
ax = fig.gca(projection='3d')

#####################
# plt.show()
# plt.savefig("SO/"+i,dpi=300)
re_bundle=re.compile(
    r"([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\d]+)[\s]+([\d]+)"
)
# url=["cluster_10up.txt","cluster_100up.txt","cluster_50up.txt","cluster_200up.txt"]
# url=["static_cluster_10up.txt","static_cluster_3up.txt"]
# url=["static_cluster_3up.txt","static_cluster_10up.txt","static_cluster_50up.txt","static_cluster_100up.txt","static_cluster_200up.txt"]
url=["static_cluster_200up.txt"]
X = []
Y = []
Z = []

X1 = []
Y1 = []
Z1 = []
X2 = []
Y2 = []
Z2 = []
X3 = []
Y3 = []
Z3 = []
X4 = []
Y4 = []
Z4 = []
X5 = []
Y5 = []
Z5 = []
X6 = []
Y6 = []
Z6 = []
X7 = []
Y7 = []
Z7 = []

R = []
G = []
B = []

# SOMA_X = []
# SOMA_Y = []
# SOMA_Z = []
QQ=0
for kk in url:
    x1 = []
    y1 = []
    z1 = []
    x2 = []
    x3 = []
    x4 = []
    y2 = []
    y3 = []
    y4 = []
    z2 = []
    z3 = []
    z4 = []
    x5 = []
    y5 = []
    z5 = []
    x6 = []
    y6 = []
    z6 = []
    x7 = []
    y7 = []
    z7 = []

    x = []
    y = []
    z = []
    cluster_list = []
    neuron_collection = []
    with open ("bundle_cluster_0p4/"+kk,"rt") as fp:
        for j in fp:
            if j[-1]=='\n':
                cluster_list.append(j[0:-1])
            else:
                cluster_list.append(j)
    for i in cluster_list:
        neuron_list = []
        with open("bundle_cluster_0p4/" + "sorting_" + i + ".txt", "rt") as fp:
            for j in fp:
                # print j
                if j[-1] == '\n':
                    j = j[:-1]
                    # k=j.split(",")
                    name = j.split("_")
                    neuron_list.append(str(nn.namelist.index(name[0])) + "_" + name[1])
                else:
                    # k=j.split(",")
                    name = j.split("_")
                    neuron_list.append(str(nn.namelist.index(name[0])) + "_" + name[1])
        # print neuron_list
        neuron_collection.append(neuron_list)
    color_index = 1
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')

    for i in neuron_collection:
        # soma_x=[]
        # soma_y=[]
        # soma_z=[]
        j = 0
        for neuron in i:
            try:
                with open('bundles_c/' + neuron + ".txt") as f:  ##############原始資料的目錄
                    j = 0
                    for line in f:  # print line
                        for iter2 in re_bundle.finditer(line):
                            # continue
                            if float(iter2.group(1))<0:
                                continue
                            else:
                                x.append(float(iter2.group(1)))
                                y.append(float(iter2.group(2)))
                                z.append(float(iter2.group(3)))
                fullneuron=neuron.split("_")
                neuron_index=int(fullneuron[0])
                # print fullneuron
                # print neuron_index
                # print nn.namelist[neuron_index]
                with open("SO_inser/SO_"+nn.namelist[neuron_index]) as fff:
                    for line in fff:
                        # print line
                        for iter2 in regex.finditer(line):
                            if float(iter2.group(1))<0:
                                continue
                            # print iter2.group(2)
                            else:
                                SO = (int(iter2.group(7)))
                                if SO == 1:
                                    x1.append(-1 * float(iter2.group(1)))
                                    y1.append( float(iter2.group(2)))
                                    z1.append( float(iter2.group(3)))
                                elif SO == 2:
                                    x2.append(-1 * float(iter2.group(1)))
                                    y2.append( float(iter2.group(2)))
                                    z2.append( float(iter2.group(3)))
                                elif SO == 3:
                                    x3.append(-1 * float(iter2.group(1)))
                                    y3.append( float(iter2.group(2)))
                                    z3.append( float(iter2.group(3)))
                                elif SO == 4:
                                    x4.append(-1 * float(iter2.group(1)))
                                    y4.append( float(iter2.group(2)))
                                    z4.append( float(iter2.group(3)))
                                elif SO == 5:
                                    x5.append(-1 * float(iter2.group(1)))
                                    y5.append( float(iter2.group(2)))
                                    z5.append( float(iter2.group(3)))
                                elif SO == 6:
                                    x6.append(-1 * float(iter2.group(1)))
                                    y6.append(float(iter2.group(2)))
                                    z6.append(float(iter2.group(3)))
                                elif SO == 7:
                                    x7.append(-1 * float(iter2.group(1)))
                                    y7.append( float(iter2.group(2)))
                                    z7.append( float(iter2.group(3)))

            except:
                print neuron
                continue
    if QQ==0:
        R.append(0.3)
        G.append(0)
        B.append(0)
    elif QQ==1:
        R.append(0.6)
        G.append(0.3)
        B.append(0)
    elif QQ==2:
        R.append(0.9)
        G.append(0.6)
        B.append(0.3)
    elif QQ==3:
        R.append(1)
        G.append(0.9)
        B.append(0.6)
    elif QQ==4:
        R.append(1)
        G.append(1)
        B.append(1)
    X.append(x)
    Y.append(y)
    Z.append(z)
    X1.append(x1)
    Y1.append(y1)
    Z1.append(z1)
    X2.append(x2)
    Y2.append(y2)
    Z2.append(z2)
    X3.append(x3)
    Y3.append(y3)
    Z3.append(z3)
    X4.append(x4)
    Y4.append(y4)
    Z4.append(z4)
    X5.append(x5)
    Y5.append(y5)
    Z5.append(z5)
    X6.append(x6)
    Y6.append(y6)
    Z6.append(z6)
    X7.append(x7)
    Y7.append(y7)
    Z7.append(z7)

    for i in xrange(len(X)):
        # plt.plot(SOMA_X[i], SOMA_Y[i], SOMA_Z[i], "o", markersize=10,color=(R[i]/10,G[i],B[i]/10,alpha))
        plt.plot(X1[i], Y1[i], Z1[i], "o", markersize=1, color=(0, 0, 1, 0.1))
        plt.plot(X2[i], Y2[i], Z2[i], "o", markersize=1, color=(0, 0.5, 0.6, 0.2))
        plt.plot(X3[i], Y3[i], Z3[i], "o", markersize=1, color=(0, 0.7, 0.4, 0.3))
        plt.plot(X4[i], Y4[i], Z4[i], "o", markersize=1, color=(0, 0.7, 0, 0.4))
        plt.plot(X5[i], Y5[i], Z5[i], "o", markersize=1, color=(0.7, 0.7, 0, 0.5))
        plt.plot(X6[i], Y6[i], Z6[i], "o", markersize=1, color=(0.8, 0.2, 0, 0.6))
        plt.plot(X7[i], Y7[i], Z7[i], "o", markersize=1, color=(1, 0, 0, 0.7))
        # plt.plot(X[i], Y[i], Z[i], "o", markersize=1, color=(0.1, 0.1, 0.1, 0.1))

    QQ=QQ+1

# x1 = []
# y1 = []
# z1 = []
# regex = re.compile(r"([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)")
# with open('soma_location_original.txt') as f:
#     for line in f:
#             # print line
#         for iter2 in regex.finditer(line):
#             x1.append(float(iter2.group(1))*10)
#             y1.append(float(iter2.group(2))*10)
#             z1.append(float(iter2.group(3))*10)
#
#
# plt.plot(x1, y1, z1, 'bo', Markersize=3)


for iii in xrange(10):
    for jjj in xrange(10):
        ax.view_init(iii * 36, jjj * 36)
        plt.savefig("bundle_cluster_0p4_fig/ALL_SO_200up_" + str(iii) + "_" + str(jjj) + ".png", dpi=100)
plt.close()
