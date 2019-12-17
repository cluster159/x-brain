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
with open("bundle_cluster_0p4/bundles_in_a_neuron.txt","wt") as fp:
    for neuron in Neuron:
        print neuron, Neuron[neuron]
        fp.writelines(nn.namelist[int(neuron)]+" "+Neuron[neuron]+"\n")
