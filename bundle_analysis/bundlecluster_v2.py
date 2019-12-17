#! /usr/bin/env python
# coding: utf-8
import numpy as np
import re
import time
import whole_brain_neuron_name as nn  # convert name to index
import os
import random
from scipy.cluster.hierarchy import linkage, dendrogram
from matplotlib import pyplot as plt
import string
import math
import pylab
import bound_list_sig20um_2sig as bound  # sphere bound
import hierachical_heatmap as heatmap
import pandas as pd

minp = bound.minp
# dendata_collection=np.matrix[20][20]

for dirPath, dirNames, fileNames in os.walk("bundle_result_long_100/"):
    for f in fileNames:
        continue
url = fileNames

x=0
y=0
z=0
cluster={}


def get_data(neuron, Time, Value):
    name=neuron.split("_")
    with open("bundle_result_long_100/" +"result_"+ str(nn.namelist.index(name[0]))+"_"+name[1]+".txt", "rt") as f_neuron:
        for line in f_neuron:
            if line.find(",")==-1:
                print line, neuron
                break
            # print line
            grouping = line.split(",")
            neuron_name = grouping[0]
            bundle_id = grouping[1][1:]
            neuron_name = neuron_name+"_"+bundle_id
            sp = grouping[2]
            vp = grouping[3][1:]
            fp = grouping[4][1:-1]
            # print neuron_name, bundle_id, sp, vp, fp, "Hi"
            # print "fp"
            # print fp
            if "e" in sp:
                sp = sp.split("e")
                sp = float(float(sp[0]) * (10 ** int(sp[1])))
            if "e" in vp:
                vp = vp.split("e")
                vp = float(float(vp[0]) * (10 ** int(vp[1])))
            if "e" in fp:
                fp = fp.split("e")
                fp = float(float(fp[0]) * (10 ** int(fp[1])))
            # print fp
            spp=int(float(sp)*100000)
            vpp=int(float(vp)*100000)
            p = int(float(fp) * 100000)
            #-------------------------
            # 演算法邏輯
            # 1. 先確認是否data就是自己，以及p是否大於閾值，是的話就不用看
            # 2. 確認data pair是否在字典裏面
            #     2-1. 若不在字典裏面→檢查相反順序的pair是否在字典裏面
            #         2-1-1. 若在的話，取min值保留，計次+1，將相似的cluster的編號全部換成最小的data的編號
            #         2-1-2. 不在的話字典加入，但不加入cluster
            #     2-2. 若在字典裡面，取min值，計次+1，將相似的cluster編號全部換成最小的data編號
             #---------------------------
            if (spp >= (0.4*100000) and vpp >= (0.7*100000) and neuron!=neuron_name):
                if (neuron,neuron_name) not in Value:
                    if (neuron_name,neuron) in Value:
                        if p < Value[(neuron_name,neuron)]:
                            Value[(neuron_name,neuron)]=p
                        Time[(neuron_name,neuron)]=Time[(neuron_name,neuron)]+1
                        # print neuron_name, neuron, Time[(neuron_name, neuron)]
                        if cluster[neuron]>cluster[neuron_name]:
                            cluster[neuron]=cluster[neuron_name]
                        else :
                            instead=cluster[neuron_name]
                            cluster[neuron_name]=cluster[neuron]
                            for n in cluster:
                                if cluster[n]==instead:
                                    cluster[n]=cluster[neuron]
                    else:
                        Value[(neuron,neuron_name)] = p
                        Time[(neuron,neuron_name)] =1
                else:
                    if p < Value[(neuron, neuron_name)]:
                        Value[(neuron, neuron_name)] = p
                    Time[(neuron, neuron_name)] = Time[(neuron, neuron_name)] + 1
                    # print neuron, neuron_name, Time[(neuron, neuron_name)]
                    if cluster[neuron] > cluster[neuron_name]:
                        cluster[neuron] = cluster[neuron_name]
                    else:
                        instead = cluster[neuron_name]
                        cluster[neuron_name] = cluster[neuron]
                        for n in cluster:
                            if cluster[n] == instead:
                                cluster[n] = cluster[neuron]
            elif p<(0.4*100000):
                break

i=0
Time = {}
Value = {}
Dist = {}
i=0
#預先給所有neuron相對應的cluster編號
#url內的是index的格式
for neuron in url:
    grouping=neuron.split("_")
    #這裡把neuron_index改成名字了
    neuron=nn.namelist[int(grouping[1])]+"_"+grouping[2][:-4]
    # print neuron
    # if i >=3:
    #     break

    # print "value"
    # print Value
    cluster[neuron]=i
    i = i + 1

i=0

#讀data
for neuron in url:
    grouping=neuron.split("_")
    neuron=nn.namelist[int(grouping[1])]+"_"+grouping[2][:-4]
    # print neuron
    # if i >=3:
    #     break
    # print i
    # print "value"
    # print Value
    get_data(neuron, Time, Value)
    i = i + 1

# print cluster


#計算每個cluster個數
count_list=[0]*(len(url))
for n in cluster:
    count_list[cluster[n]]=count_list[cluster[n]]+1



##寫出>#的cluster有哪些
for j in xrange(len(count_list)):
    if count_list[j] >= 200:
        with open("bundle_cluster_sp4_vp7/static_cluster_200up.txt", "at") as fp:
            fp.writelines("cluster_" + str(j) + '\n')
    elif count_list[j] >= 100:
        with open("bundle_cluster_sp4_vp7/static_cluster_100up.txt", "at") as fp:
            fp.writelines("cluster_" + str(j) + '\n')
    elif count_list[j] >= 50:
        with open("bundle_cluster_sp4_vp7/static_cluster_50up.txt", "at") as fp:
            fp.writelines("cluster_" + str(j) + '\n')
    elif count_list[j] >= 10:
        with open("bundle_cluster_sp4_vp7/static_cluster_10up.txt", "at") as fp:
            fp.writelines("cluster_" + str(j) + '\n')
    elif count_list[j] >= 3:
        with open("bundle_cluster_sp4_vp7/static_cluster_3up.txt", "at") as fp:
            fp.writelines("cluster_" + str(j) + '\n')
    elif count_list[j] >= 1:
        with open("bundle_cluster_sp4_vp7/static_cluster_1up.txt", "at") as fp:
            fp.writelines("cluster_" + str(j) + '\n')

# ##寫出每個cluster的成員--
# for n in cluster:
#     with open ("bundle_cluster_sp4_vp7/cluster_"+str(cluster[n])+".txt","at") as fp:
#         grouping=n.split("_")
#         neuron_index=nn.namelist.index(grouping[0])
#         fp.writelines(str(grouping[0])+", "+str(neuron_index)+"_" +str(grouping[1])+ "\n")

##
assembly_collection=[]
assembly_index=[]
for i in xrange(len(count_list)):
    if count_list[i]  > 2 :
        assembly=[]
        for j in cluster:
            if cluster[j]==i:
                assembly.append(j)
        assembly_collection.append(assembly)
        assembly_index.append(i)
print assembly_collection

for j in xrange(len(assembly_collection)):
    for w in assembly_collection[j]:
        name = w.split("_")
        # print "name:",name
        with open("bundle_result_long_100/" + "result_" + str(nn.namelist.index(name[0])) + "_" + name[
            1] + ".txt", "rt") as f_neuron:
            for line in f_neuron:
                if line.find(",") == -1:
                    print line, neuron
                    break
                # print line
                grouping = line.split(",")
                neuron_name = grouping[0]
                bundle_id = grouping[1][1:]
                neuron_name = neuron_name + "_" + bundle_id
                sp = grouping[2]
                vp = grouping[3][1:]
                fp = grouping[4][1:-1]
                # print neuron_name, bundle_id, sp, vp, fp, "Hi"
                # print "fp"
                # print fp
                if "e" in sp:
                    sp = sp.split("e")
                    sp = float(float(sp[0]) * (10 ** int(sp[1])))
                if "e" in vp:
                    vp = vp.split("e")
                    vp = float(float(vp[0]) * (10 ** int(vp[1])))
                if "e" in fp:
                    fp = fp.split("e")
                    fp = float(float(fp[0]) * (10 ** int(fp[1])))
                # print fp
                p = int(float(fp)*100000)
                if neuron_name in assembly_collection[j]:
                    if (neuron_name, w) in Value:
                        if Time[(neuron_name, w)] > 1:
                            continue
                        else:
                            if p == Value[(neuron_name, w)]:
                                continue
                            else:
                                # print "one->two_1", p
                                if p < Value[(neuron_name, w)]:
                                    Value[(neuron_name, w)] = p
                                    Time[(neuron_name, w)] = Time[(neuron_name, w)] + 1
                                else:
                                    Time[(neuron_name, w)] = Time[(neuron_name, w)] + 1
                    elif (w, neuron_name) in Value:
                        if Time[(w, neuron_name)] > 1:
                            continue
                        else:
                            # print "one->two_2", p
                            if p == Value[(w, neuron_name)]:
                                continue
                            else:
                                if p < Value[(w, neuron_name)]:
                                    Value[(w, neuron_name)] = p
                                    Time[(w, neuron_name)] = Time[(w, neuron_name)] + 1
                                else:
                                    Time[(w, neuron_name)] = Time[(w, neuron_name)] + 1
                    elif (neuron_name, w) not in Value:
                        Value[(neuron_name, w)] = p
                        Time[(neuron_name, w)] = 1
                        if neuron_name==w:
                            Time[(neuron_name, w)]=2
                        # print "both_not", p

##找出每個cluster成員，任兩個的關係

for j in xrange(len(assembly_collection)):
    with open ("bundle_cluster_sp4_vp7/cluster_"+str(assembly_index[j])+"_correlation_map.txt","wt") as fp:
        fp.writelines("cluster_"+str(assembly_index[j]) + " ")
        line=""
        for w in assembly_collection[j]:
            line=line+" "+w
        fp.writelines(line+"\n")
        for w in assembly_collection[j]:
            fp.writelines(str(w)+" ")
            line=""
            for m in assembly_collection[j]:
                if (m, w) in Value:
                    if Time[(m, w)] > 1:
                        line = line + str(float(Value[(m, w)])/100000) + " "
                    else:
                        line = line + "0 "
                elif (w, m) in Value:
                    if Time[(w, m)] > 1:
                        line = line + str(float(Value[(w, m)])/100000) + " "
                    else:
                        line = line + "0 "
                elif m == w:
                    line = line + "1 "
                else:
                    line = line + "-1 "
                if (m, w) in Value and (w, m) in Value and m!=w:
                    print "error\n"
                # if (m, w) not in Value and (w, m) not in Value:
                #     print "Notice\n"
            fp.writelines(line + "\n")



