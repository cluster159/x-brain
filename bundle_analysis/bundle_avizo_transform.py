#! /usr/bin/env python
# coding: utf-8


import numpy as np
# import pandas as pd
import re
# from itertools import cycle
# import  operator
from collections import OrderedDict
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pylab

import cPickle as pickle  # needed while saving data
import time
import whole_brain_neuron_name as nn
import os
import random

start = time.time()
xyz = {}


class Scalar:
    value = 0.5
    sigma = 3.0


def radius(x, y, z, xm, ym, zm):  # distance between two points
    r = np.sqrt((x - xm) ** 2 + (y - ym) ** 2 + (z - zm) ** 2)
    return r


def prob(r, sigma):  # probability formula
    p = np.exp(-(r / sigma) ** 2 / 2) / (np.sqrt(2 * (np.pi) ** 3) * (sigma) ** 3)
    return p


def toScalar(x):  # convert to resolution
    if round((x % Scalar.value) / Scalar.value) == 1:
        x = x - x % Scalar.value + Scalar.value
    else:
        x = x - x % Scalar.value
    return float(x)  # actually float(x) already present the same value as float("%.3f" %x)


def name_to_index(name):  # convet name to int
    return nn.namelist.index(name)


def pos_to_int(x):  # convert to int(x*10)
    return int(round(10 * x))


def theta_to_int(x):  # convert to int(x*100)
    return int(round(100 * x))


countpoints = 0
allpoints = OrderedDict()
exist = set()  # ***check if exist is needed to be outside the for loop
xyz = []
time = 0
swc = " "
count = 0
# ############somata
x1 = []
y1 = []
z1 = []
regex = re.compile(r"([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)")
with open('soma_location_original.txt') as f:
    for line in f:
        # print line
        for iter2 in regex.finditer(line):
            x1.append(float(iter2.group(1)) * 10)
            y1.append(float(iter2.group(2)) * 10)
            z1.append(float(iter2.group(3)) * 10)

fig = plt.figure()
ax = fig.gca(projection='3d')

plt.plot(x1, y1, z1, 'bo', Markersize=3)
#####################
# plt.show()
# plt.savefig("SO/"+i,dpi=300)
re_bundle = re.compile(
    r"([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\d]+)[\s]+([\d]+)"
)
# url=["cluster_10up.txt","cluster_100up.txt","cluster_50up.txt","cluster_200up.txt"]
# url=["static_cluster_10up.txt","static_cluster_3up.txt"]
url = ["static_cluster_3up.txt", "static_cluster_10up.txt", "static_cluster_50up.txt", "static_cluster_100up.txt",
       "static_cluster_200up.txt"]
path_record="FC_bundle_Avizo_0p4_all/"
if not os.path.isdir(path_record):
    os.mkdir(path_record)


for kk in url:
    # print kk
    x = []
    y = []
    z = []
    cluster_list = []
    neuron_collection = []
    with open("bundle_cluster_0p4/" + kk, "rt") as fp:
        for j in fp:
            if j[-1] == '\n':
                cluster_list.append(j[0:-1])
            else:
                cluster_list.append(j)
#####################
    xyz=""
    length=0
#####################
    for i in cluster_list:
        neuron_list = []
        with open("bundle_cluster_0p4/" + "sorting_" + i + ".txt", "rt") as fp:
            for j in fp:
                # print j
                if j[-1] == '\n':
                    j = j[:-1]
                    name = j.split("_")
                    neuron_list.append(str(nn.namelist.index(name[0])) + "_" + name[1])
                else:
                    name = j.split("_")
                    neuron_list.append(str(nn.namelist.index(name[0])) + "_" + name[1])
        # print i
#############################
        # length = 0
        # xyz = ""
#############################
        for neuron in neuron_list:
            # print neuron
            try:
                with open('bundles_c/' + neuron + ".txt") as f:  ##############原始資料的目錄
                    j = 0
                    for line in f:  # print line
                        for iter2 in re_bundle.finditer(line):
                            xyz = xyz + str(float(iter2.group(1)) / 10) + " " + str(
                                float(iter2.group(2)) / 10) + " " + str(float(iter2.group(3)) / 10) + "\n"
                            length = length + 1
            except:
                print "ERROR"
                # print neuron
                break
###################################
        # record_lines = ""
        # record_lines = "# AmiraMesh 3D ASCII 2.0\n" + "define Markers " + str(
        #     length) + "\nParameters {\n" + "    ContentType \"LandmarkSet\",\n" + "    NumSets 1\n" + "}\n" + "Markers { float[3] Coordinates } @1\n" + "# Data section follows\n" + "@1\n" + xyz
        # with open(path_record + kk[15:-4] + "_" + i + ".am", "wt")as f:
        #     f.write(record_lines)
####################################
    record_lines = ""
    record_lines = "# AmiraMesh 3D ASCII 2.0\n" + "define Markers " + str(
        length) + "\nParameters {\n" + "    ContentType \"LandmarkSet\",\n" + "    NumSets 1\n" + "}\n" + "Markers { float[3] Coordinates } @1\n" + "# Data section follows\n" + "@1\n" + xyz
    with open(path_record + kk[15:-4] + ".am", "wt")as f:
        f.write(record_lines)