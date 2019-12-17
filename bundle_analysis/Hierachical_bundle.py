#! /usr/bin/env python
# coding: utf-8
import numpy as np
import re
import time
import whole_brain_neuron_name as nn #convert name to index
import os
import random
from scipy.cluster.hierarchy import linkage, dendrogram
from matplotlib import pyplot as plt
import string
import math
import pylab
import bound_list_sig20um_2sig as bound #sphere bound
from scipy.cluster.hierarchy import cophenet
from scipy.cluster.hierarchy import linkage, dendrogram,to_tree,leaves_list,leaders,fcluster
from scipy.spatial.distance import pdist
import pandas as pd
import seaborn as sns
sns.set()
minp=bound.minp
# dendata_collection=np.matrix[20][20]
url_list=[]
pathway="bundle_cluster_0p4/"
for dirPath, dirNames, fileNames in os.walk(pathway):
    for f in fileNames:
        if f.find("up")!=-1 and f.find("1up")==-1:
            url_list.append(f)
cluster_list = []

for list in url_list:
    with open (pathway+list,"rt") as f:
        for line in f:
            cluster_list.append(line[:-1])
    # print cluster_list
# print cluster_list_collection
subfix="_correlation_map.txt"


for cluster in cluster_list:
    Distance_matrix=[]
    Heat_matrix=[]
    name=[]
    with open(pathway+cluster+subfix,"rt") as f:
        i=0
        for line in f:
            j=0
            grouping=line.split(" ")

            # print grouping
            distance_row = []
            heat_row=[]
            if i==0:
                leaf_number=len(grouping)-2
                # print "number", leaf_number
                i=1
                grouping[leaf_number+1]=grouping[leaf_number+1][:-1]
                name=grouping[2:]
                continue
            else:
                for m in grouping[1:-1]:
                    if float(m)==0 or float(m)==-1:
                        distance=100000
                        heat=0
                    else:
                        distance=1/float(m)
                        heat=float(m)
                    distance_row.append(distance)
                    heat_row.append(heat)
            # print distance_row
            Distance_matrix.append(distance_row)
            Heat_matrix.append(heat_row)
    Distance_matrix=np.array(Distance_matrix)
    linkage_matrix = linkage(Distance_matrix, "single")
    new_ids = leaves_list(linkage_matrix)
    new_name=[]
    new_name_string=""
    for i in xrange(leaf_number):
        new_name.append(name[new_ids[i]])
        new_name_string=new_name_string+name[new_ids[i]]+"\n"
    # print new_name
    # print name
    # print new_ids
    with open(pathway+"sorting_"+cluster+".txt","wt") as f:
        f.writelines(new_name_string)

    dendrogram(linkage_matrix)
    # plt.show()
    plt.savefig(pathway[:-1]+"_fig/"+cluster+"_single.png",dpi=100)
    plt.close()
    #######Sorting###########
    Heat_new_matrix=[]
    for i in new_ids:
        heat_new_row = []
        for j in new_ids:
            heat_new_row.append(Heat_matrix[i][j])
        Heat_new_matrix.append(heat_new_row)
    Heat_new_matrix=np.array(Heat_new_matrix)
    ###########################
    fig, ax = plt.subplots()
    im = ax.imshow(Heat_new_matrix)
    cbar = ax.figure.colorbar(im,label="Similarity",ticks=[0.0,0.2,0.4,0.6,0.8,1.0])
    # ax.set_clim([0.0,1.0])
    # plt.show()
    plt.savefig(pathway[:-1]+"_fig/"+cluster+"_heatmap_single.png",dpi=100)
    plt.close()
