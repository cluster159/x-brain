#! /usr/bin/env python
# coding: utf-8

#####z覆蓋到了..............................................................................
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

countpoints = 0
allpoints = OrderedDict()
exist = set()  # ***check if exist is needed to be outside the for loop
xyz=[]
time = 0
swc = " "
count=0

x1 = []
y1 = []
z1 = []
regex = re.compile(r"([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)")
re_bundle=re.compile(r"([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\d]+)")
# url=["cluster_10up.txt","cluster_100up.txt","cluster_50up.txt","cluster_200up.txt"]
url=["static_cluster_200up.txt"]
# url=[]
cluster_list=[]
neuron_collection=[]
for i in url:
    with open ("xbundle_22500_sig30_50_cluster_0p4/"+i,"rt") as fp:
        for j in fp:
            if j[-1]=='\n':
                cluster_list.append(j[0:-1])
            else:
                cluster_list.append(j)
for i in cluster_list:
    neuron_list=[]
    try:
        with open("xbundle_22500_sig30_50_cluster_0p4/sorting_" + i + ".txt", "rt") as fp:
            for j in fp:
                if j[-1] == '\n':
                    j = j[:-1]
                    neuron_list.append(j)
                else:
                    neuron_list.append(j)
        neuron_collection.append(neuron_list)
    except:
        print "error", i
        continue
R=[]
G=[]
B=[]
alpha=0.7
X=[]
Y=[]
Z=[]
# SOMA_X=[]
# SOMA_Y=[]
# SOMA_Z=[]
# fig = plt.figure()
# ax = fig.gca(projection='3d')
#
for i in neuron_collection:
    r=random.random()
    g=random.random()
    b=random.random()
    x=[]
    y=[]
    z=[]
    # soma_x=[]
    # soma_y=[]
    # soma_z=[]
    j = 0
    for neuron in i:
        with open('x_bundles_c_0p3_linux/' + neuron + ".txt") as f:  ##############原始資料的目錄
            j = 0
            for line in f:  # print line
                for iter2 in re_bundle.finditer(line):
                    x.append(int(iter2.group(1)))
                    y.append(int(iter2.group(2)))
                    z.append(int(iter2.group(3)))
                    # if j == 0:
                    #     soma_x.append(float(iter2.group(1)))
                    #     soma_y.append(float(iter2.group(2)))
                    #     soma_z.append(float(iter2.group(3)))
                    #     j = j + 1
        # plt.plot(x, y, z, "o", markersize=3, color=(random.random(), random.random(), random.random(), alpha))
    # plt.show()
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')

    X.append(x)
    Y.append(y)
    Z.append(z)
    # SOMA_X.append(soma_x)
    # SOMA_Y.append(soma_y)
    # SOMA_Z.append(soma_z)
    R.append(r)
    G.append(g)
    B.append(b)

fig = plt.figure()
ax = fig.gca(projection='3d')

for i in xrange(len(X)):
    # print R[i],G[i],B[i]
    # plt.plot(SOMA_X[i], SOMA_Y[i], SOMA_Z[i], "o", markersize=10,color=(R[i]/10,G[i],B[i]/10,alpha))
    plt.plot(X[i], Y[i], Z[i], "o", markersize=6,color=(R[i],G[i],B[i],alpha))
# ------------------------------------------------------
# 20170727_24_3-brain_charng.swc_0_0
try:
    os.mkdir("xbrain_22500_bundlefig1")
except:
    pass
re_swc = re.compile(r"([\d]+)[\s]+([\d]+)[\s]+([\d]+[\.][\d]+)[\s]+([\d]+[\.][\d]+)[\s]+([\d]+[\.][\d]+)[\s]+([\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+)")
path_original="x_brain/22500_swc_final/"
x_original = []
y_original = []
z_original = []
with open(path_original + "20170727_24_3-brain_charng.swc", "rt")as f:
    for line in f:
        for iter in re_swc.finditer(line):
            x_original.append(float(iter.group(3))*3)
            y_original.append(float(iter.group(4))*3)
            z_original.append(float(iter.group(5))*3)
plt.plot(x_original, y_original, z_original, "ko", markersize=2,alpha=0.3)
for iii in xrange(6):
    for jjj in xrange(6):
        ax.view_init(iii * 60, jjj * 60)
        fig.set_size_inches(18.5, 18.5)
        plt.savefig("xbrain_22500_bundlefig1/sig30_200up_" + str(iii) + "_" + str(jjj) + ".png",dpi=100)
plt.close()
# plt.show()
