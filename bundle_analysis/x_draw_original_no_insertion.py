#! /usr/bin/env python
# coding: utf-8
import numpy as np
import re
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pylab
import cPickle as pickle #needed while saving data
import time
import os
import random
start = time.time()
xyz = {}
class Scalar:
    value = 0.5
    sigma = 3.0


path_original="x_brain/22500_swc_final/"
path_without_insertion="x_bundles_c_no_insertion_0p3_linux/"
re_swc = re.compile(r"([\d]+)[\s]+([\d]+)[\s]+([\d]+[\.][\d]+)[\s]+([\d]+[\.][\d]+)[\s]+([\d]+[\.][\d]+)[\s]+([\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+)")
re_without_insertion=re.compile(r"([\d]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\d]+)")

warping_list=[]
with open("x_swc_warped_22500_list.txt","rt") as f:
    for line in f:
        if line[-1]=="\n":
            line=line[:-1]
        warping_list.append(line)

prefix=[]
suffix=[]
os.mkdir("xbrain_22500_brainfig1")

# for dirPath, dirNames, fileNames in os.walk(path_without_insertion):
#     continue
# print fileNames[0]
#
# for name in fileNames:
#     grouping = name.split("_")
#     prefix.append(int(grouping[0]))
#     suffix.append(grouping[1])
# print prefix[0]
for brain_index in xrange(len(warping_list)):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    x_original = []
    y_original = []
    z_original = []
    with open (path_original+warping_list[brain_index],"rt")as f:
        for line in f:
            for iter in re_swc.finditer(line):
                x_original.append(float(iter.group(3)))
                y_original.append(float(iter.group(4)))
                z_original.append(float(iter.group(5)))
    plt.plot(x_original,y_original,z_original,"ko",markersize=5)
    # check=0
    # for i in xrange(len(prefix)):
        # x_fragmentize = []
        # y_fragmentize = []
        # z_fragmentize = []
        # if prefix[i]==brain_index:
        #     with open(path_without_insertion+str(prefix[i])+"_"+suffix[i],"rt") as f:
        #         for line in f:
        #             for iter in re_without_insertion.finditer(line):
        #                 x_fragmentize.append(float(iter.group(1))/3)
        #                 y_fragmentize.append(float(iter.group(2))/3)
        #                 z_fragmentize.append(float(iter.group(3))/3)
        #     check=1
        # if check ==1 and prefix[i]!=brain_index:
        #     print prefix[i], brain_index
        #     break
        # plt.plot(x_fragmentize,y_fragmentize,z_fragmentize,"r-",linewidth="3")
    # plt.show()
    for iii in xrange(6):
        for jjj in xrange(6):
            ax.view_init(iii * 60, jjj * 60)
            fig.set_size_inches(18.5, 18.5)
            plt.savefig("xbrain_22500_brainfig1/"+warping_list[brain_index] + "_" + str(iii) + "_" + str(jjj) + ".png", dpi=100)
    plt.close()