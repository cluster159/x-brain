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

regex = re.compile(
    r"([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\d]+)"
)

countpoints = 0
allpoints = OrderedDict()
exist = set()  # ***check if exist is needed to be outside the for loop
xyz=[]
time = 0
swc = " "
count=0
re_bundle=re.compile(
    r"([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\d]+)[\s]+([\d]+)"
)
# url=["cluster_3up.txt","cluster_10up.txt","cluster_50up.txt","cluster_100up.txt","cluster_200up.txt"]
url=["static_cluster_3up.txt","static_cluster_10up.txt","static_cluster_50up.txt","static_cluster_100up.txt","static_cluster_200up.txt"]

Neuron={}
cluster_collection=[]
neuron_collection=[]
for i in url:
    cluster_list = []
    with open ("bundle_cluster_0p4/"+i,"rt") as fp:
        for j in fp:
            if j[-1]=='\n':
                cluster_list.append(j[0:-1])
            else:
                cluster_list.append(j)
    cluster_collection.append(cluster_list)
list=[3,10,50,100,200]
list_index=0
for cluster_list in cluster_collection:
    for i in cluster_list:
        neuron_list = []
        with open("bundle_cluster_0p4/" + "sorting_" + i + ".txt", "rt") as fp:
            for j in fp:
                # print j
                if j[-1] == '\n':
                    j = j[:-1]
                    # k=j.split(",")
                    name = j.split("_")
                    if str(nn.namelist.index(name[0])) not in Neuron:
                        Neuron[str(nn.namelist.index(name[0]))]=i+"_"+name[1]+"_"+str(list[list_index])
                    #     cluster+branch_id+number
                    else:
                        Neuron[str(nn.namelist.index(name[0]))] =Neuron[str(nn.namelist.index(name[0]))]+","+i+"_"+name[1]+"_"+str(list[list_index])
                    neuron_list.append(str(nn.namelist.index(name[0])) + "_" + name[1])
                else:
                    # k=j.split(",")
                    name = j.split("_")
                    if str(nn.namelist.index(name[0])) not in Neuron:
                        Neuron[str(nn.namelist.index(name[0]))]=i+"_"+name[1]+"_"+str(list[list_index])
                    #     cluster+branch_id+number
                    else:
                        Neuron[str(nn.namelist.index(name[0]))] =Neuron[str(nn.namelist.index(name[0]))]+","+i+"_"+name[1]+"_"+str(list[list_index])

                    neuron_list.append(str(nn.namelist.index(name[0])) + "_" + name[1])
        # print neuron_list
        neuron_collection.append(neuron_list)
    list_index=list_index+1

###############################draw
R=[]
G=[]
B=[]
alpha=0.7

# Xn=[]
# Yn=[]
# Zn=[]

X=[]
Y=[]
Z=[]
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

SOMA_X=[]
SOMA_Y=[]
SOMA_Z=[]
color_index=1
fig = plt.figure()
ax = fig.gca(projection='3d')
for neuron in Neuron:
    cluster_number = Neuron[neuron].count(",") + 1
    if Neuron[neuron].find("_188_")==-1:
        continue
    # print "HERE"
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
    soma_x = []
    soma_y = []
    soma_z = []
    r = random.random()
    g = random.random()
    b = random.random()
    # print nn.namelist[int(neuron)]
    soma_check = 0
    with open("SO_inser/SO_" + nn.namelist[int(neuron)]) as fff:
        # print "SUCCESS1"
        for line in fff:
            # print line
            for iter2 in regex.finditer(line):
                if soma_check < 2:
                    soma_x.append(float(iter2.group(1)))
                    soma_y.append(float(iter2.group(2)))
                    soma_z.append(float(iter2.group(3)))
                    soma_check = soma_check + 1

                SO = (int(iter2.group(7)))
                if SO == 1:
                    x1.append(float(iter2.group(1)))
                    y1.append(float(iter2.group(2)))
                    z1.append(float(iter2.group(3)))
                elif SO == 2:
                    x2.append(float(iter2.group(1)))
                    y2.append(float(iter2.group(2)))
                    z2.append(float(iter2.group(3)))
                elif SO == 3:
                    x3.append(float(iter2.group(1)))
                    y3.append(float(iter2.group(2)))
                    z3.append(float(iter2.group(3)))
                elif SO == 4:
                    x4.append(float(iter2.group(1)))
                    y4.append(float(iter2.group(2)))
                    z4.append(float(iter2.group(3)))
                elif SO == 5:
                    x5.append(float(iter2.group(1)))
                    y5.append(float(iter2.group(2)))
                    z5.append(float(iter2.group(3)))
                elif SO == 6:
                    x6.append(float(iter2.group(1)))
                    y6.append(float(iter2.group(2)))
                    z6.append(float(iter2.group(3)))
                elif SO == 7:
                    x7.append(float(iter2.group(1)))
                    y7.append(float(iter2.group(2)))
                    z7.append(float(iter2.group(3)))
    plt.plot(x1, y1, z1, "o", markersize=1, color=(r, g, b, 0.9))
    plt.plot(x2, y2, z2, "o", markersize=1, color=(r, g, b, 0.8))
    plt.plot(x3, y3, z3, "o", markersize=1, color=(r, g, b, 0.7))
    plt.plot(x4, y4, z4, "o", markersize=1, color=(r, g, b, 0.6))
    plt.plot(x5, y5, z5, "o", markersize=1, color=(r, g, b, 0.5))
    plt.plot(x6, y6, z6, "o", markersize=1, color=(r, g, b, 0.4))
    plt.plot(x7, y7, z7, "o", markersize=1, color=(r, g, b, 0.3))
    plt.plot(soma_x, soma_y, soma_z, "o", markersize=15, color=(r, g, b, 1))

    ####################draw_bundle_branch
    cluster_number = Neuron[neuron].count(",") + 1
    print Neuron[neuron]
    print cluster_number
    cluster_number = Neuron[neuron].count(",") + 1

    if cluster_number > 1:
        cluster_sequence = Neuron[neuron].split(",")
        for i in cluster_sequence:
            x = []
            y = []
            z = []
            grouping = i.split("_")
            branch_id = grouping[2]
            cluster_id = grouping[1]
            size = int(grouping[3])
            neuron_index = neuron
            print i
            print grouping
            print neuron_index + "_" + branch_id
            with open('bundles_c/' + str(neuron_index) + "_" + str(branch_id) + ".txt") as f:  ##############原始資料的目錄
                print "SUCCESS2"
                j = 0
                for line in f:  # print line
                    for iter2 in re_bundle.finditer(line):
                        x.append(float(iter2.group(1)))
                        y.append(float(iter2.group(2)))
                        z.append(float(iter2.group(3)))
            if size == 200:
                branch_size = 6
            elif size == 100:
                branch_size = 5
            elif size == 50:
                branch_size = 4
            elif size == 10:
                branch_size = 3
            elif size == 3:
                branch_size = 2
            plt.plot(x, y, z, "o", markersize=branch_size + 4, color=(r, g, b, 0.2))


    elif cluster_number==1:
        x = []
        y = []
        z = []
        grouping = Neuron[neuron].split("_")
        branch_id = grouping[2]
        cluster_id = grouping[1]
        size = int(grouping[3])
        neuron_index = neuron
        print i
        print grouping
        print neuron_index + "_" + branch_id
        with open('bundles_c/' + str(neuron_index) + "_" + str(branch_id) + ".txt") as f:  ##############原始資料的目錄
            print "SUCCESS3"
            j = 0
            for line in f:  # print line
                for iter2 in re_bundle.finditer(line):
                    x.append(float(iter2.group(1)))
                    y.append(float(iter2.group(2)))
                    z.append(float(iter2.group(3)))
        if size == 200:
            branch_size = 6
        elif size == 100:
            branch_size = 5
        elif size == 50:
            branch_size = 4
        elif size == 10:
            branch_size = 3
        elif size == 3:
            branch_size = 2
        plt.plot(x, y, z, "o", markersize=branch_size+4, color=(r, g, b, 0.2))


        # for iii in xrange(10):
        #     for jjj in xrange(10):
        #         fig.set_size_inches(18.5, 18.5)
        #         ax.view_init(iii * 36, jjj * 36)
        #         plt.savefig(
        #             "bundle_cluster_0p4_fig/NEURON_100up" + cluster_list[neuron_collection.index(i)] + "_" + str(
        #                 iii) + "_" + str(jjj) + ".png", dpi=100)

        # plt.show()
        # plt.close()
    plt.show()
    fig = plt.figure()
    ax = fig.gca(projection='3d')

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
    SOMA_X.append(soma_x)
    SOMA_Y.append(soma_y)
    SOMA_Z.append(soma_z)
    R.append(r)
    G.append(g)
    B.append(b)
plt.show()
fig = plt.figure()
ax = fig.gca(projection='3d')
# for iii in xrange(10):
#     for jjj in xrange(10):
#         fig.set_size_inches(18.5, 18.5)
#         ax.view_init(iii * 36, jjj * 36)
#         plt.savefig("bundle_cluster_0p4_fig/NEURON_100up" + cluster_list[neuron_collection.index(i)] + "_" + str(
#             iii) + "_" + str(jjj) + ".png", dpi=100)

    # plt.show()
plt.close()
