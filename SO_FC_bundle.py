#! /usr/bin/env python
# coding: utf-8
import numpy as np
import re
import time
import SO_FC_bundle_spatial as s_infer
import SO_FC_bundle_vector as v_infer
import bound_list_sig20um_2sig as bound #sphere bound
import os
import random
import math
import threading
import whole_brain_neuron_name as nn
from collections import OrderedDict

start = time.time()
resolution = 0.5
cube_size=4
path_neuron="flycircuit_branch/"
similarity_result_path="result_bundle_SO_branch_similariy0p4_sig2um/"
maxp = bound.maxp
minp = bound.minp #sigma=32取2個sigma
stdp = math.log10(maxp/minp)
def to_int(original_number):
    if original_number<0:
        mark=-1
        original_number=original_number*-1
    else:
        mark=1
    k=(original_number-int(original_number))
    if k<0.25:
        original_number=int(original_number)* 10
    elif k<0.75:
        original_number = int(original_number) * 10 + 5
    else:
        original_number = int(original_number) * 10 +10
    return original_number*mark

def pos_to_int(x): #convert to int(x*10)
    return int(round(10*x))

def name_to_index(name): #the input name should be in "filename" form
    return nn.namelist.index(name)

detect_list = []  # list: the sphere bound for the input point→建一個以input point為中心的球狀座標表單
for j in bound.bounding_list:
    detect_list.append([pos_to_int(j['x']),pos_to_int(j['y']),pos_to_int(j['z']),math.log10(j['p'] / minp) / stdp])
detect_list = sorted(detect_list, reverse=True, key=lambda x: x[3])
print detect_list[0]

neuron_list=os.listdir(path_neuron)

if not os.path.isdir(similarity_result_path):
    os.mkdir(similarity_result_path)

def add_neuro_space(points,brain_dictionary):
    # print "load neuron into background"
    # print points

    for node in points:
        # print node
        x, y, z = node[0], node[1], node[2]
        if x > 0:
            x_index = int(abs(x) / cube_size)
        else:
            x_index = -1 * int(abs(x) / cube_size)
        if y > 0:
            y_index = int(abs(y) / cube_size)
        else:
            y_index = -1 * int(abs(y) / cube_size)
        if z > 0:
            z_index = int(abs(z) / cube_size)
        else:
            z_index = -1 * int(abs(z) / cube_size)
        if (x_index, y_index, z_index) in brain_dictionary:
            if (x,y,z) not in brain_dictionary[(x_index, y_index, z_index)]:
                brain_dictionary[(x_index, y_index, z_index)][(x,y,z)]=[(node[3],node[4],node[5],node[6],node[7],node[8])]
            else:
                brain_dictionary[(x_index, y_index, z_index)][(x,y,z)].append((node[3],node[4],node[5],node[6],node[7],node[8]))
        else:
            brain_dictionary[(x_index, y_index, z_index)] = {}
            brain_dictionary[(x_index, y_index, z_index)][(x, y, z)] = [(node[3], node[4], node[5], node[6], node[7], node[8])]
    return brain_dictionary

brain_dictionary={}
trunks = []

for neuron in neuron_list:
    print "load neuron trunks", neuron
    # print "brain_dictionary"
    # print brain_dictionary
    raw_points=[]
    points=[]
    old_neuron_id = -1
    old_branch_id = -1
    start_node = -1
    with open(path_neuron + neuron) as f:
        neuron_id = nn.namelist.index(neuron)
        for line in f:
            # print line
            if line[-1]=='\n':
                line=line[:-1]
            grouping = line.split(" ")
            # x,y,z,branch_id,node_id,parent,SO
            # print grouping
            x, y, z, branch_id, node_id, parent, SO = float(grouping[0]), float(grouping[1]), float(grouping[2]), int(grouping[3]), int(grouping[4]), int(grouping[5]),int(grouping[6])
            if not branch_id==old_branch_id:
            # if not neuron_id == old_neuron_id:
                if len(points)!=0:
                    trunks.append(points)
                    brain_dictionary = add_neuro_space(points, brain_dictionary)

                    points=[]
                    raw_points=[]
                old_branch_id, old_neuron_id = branch_id, neuron_id
                raw_points.append([x, y, z, neuron_id, branch_id, node_id, parent, SO])
                # points.append([to_int(x), to_int(y), to_int(z), neuron_id, branch_id, node_id, 0, 0, 0])
                ##i如果是單條trunk 直接一直接前一個點就好~如果是多個trunk....所以可能要用搜尋的方式比較快
            else:
                raw_points.append([x, y, z, neuron_id, branch_id, node_id, parent, SO])
                tmp_order=len(raw_points)
                tmp_order=tmp_order-1
                # print tmp_order
                find=0
                count=0
                while tmp_order >= 0:
                    # print tmp_order
                    # count=count+1
                    # print count
                    if parent!=raw_points[tmp_order][5]:
                        tmp_order=tmp_order-1
                    else:
                        find=1
                        break
                if find==0:
                    print "error"
                    break
                vx = x - raw_points[tmp_order][0]
                vy = y - raw_points[tmp_order][1]
                vz = z - raw_points[tmp_order][2]
                r = (vx**2 + vy**2 +vz**2)**0.5
                if r != 0:
                    vx, vy, vz = vx/r, vy/r, vz/r
                else:
                    print "r==0"
                insertion_num = int(r / resolution) ## resolution = 0.5
                # if points[tmp_order][6]==0 and points[tmp_order][7]==0 and
                for insertion in xrange(0,insertion_num):
                    insertion_x, insertion_y, insertion_z = raw_points[tmp_order][0]+resolution*vx*insertion, raw_points[tmp_order][1]+resolution*vy*insertion, raw_points[tmp_order][2]+resolution*vz*insertion
                    if ((insertion_x-raw_points[tmp_order][0])**2+(insertion_x-raw_points[tmp_order][0])**2+(insertion_x-raw_points[tmp_order][0])**2)**0.5>=r:
                        break
                    insertion_x, insertion_y, insertion_z = to_int(insertion_x), to_int(insertion_y), to_int(insertion_z)
                    # if insertion_x == -750:
                        # print insertion_x, insertion_y, insertion_z, insertion_num
                        # os.system("pause")
                    if insertion == 0:
                        points.append([insertion_x, insertion_y, insertion_z, neuron_id, branch_id, node_id, vx, vy, vz])
                    elif not (insertion_x == points[-1][0] and insertion_y ==points[-1][1] and insertion_z == points[-1][2]):
                        points.append([insertion_x, insertion_y, insertion_z, neuron_id, branch_id, node_id, vx, vy, vz])
                # if not (insertion_x == to_int(x) and insertion_y == y and insertion_z == z):
                #     points.append([insertion_x, insertion_y, insertion_z, neuron_id, branch_id, node_id, vx, vy, vz])
                # print points
    if len(points)!=0:
        trunks.append(points)
        brain_dictionary = add_neuro_space(points,brain_dictionary)
    if find ==0:
        print "TMP"
        print nn.namelist[neuron_id]
        break

# print "trunks"
# print trunks
# print "brain"
# print brain_dictionary

# print brain_dictionary

if find==0:
    os.system("pause")

#
# ################
def similarity_calculation(brain_dictionary,trunks,detect_list, cube_size):
    full_len=float(len(trunks))
    tmp_count=0
    for inputpoints in trunks:
        print "progress:",tmp_count,full_len
        tmp_count=tmp_count+1
        final = []
        s_inference = s_infer.inference(inputpoints=inputpoints, detect_list=detect_list,
                                      brain_dictionary=brain_dictionary, cube_size=cube_size)
        # print "s",s_inference
        v_inference = v_infer.inference(inputpoints=inputpoints, detect_list=detect_list,
                                      brain_dictionary=brain_dictionary, cube_size=cube_size)
        # print "v",v_inference

        final = []
        for n in s_inference:
            # print n
            if n not in v_inference:
                v_inference[n] = 0
            final.append((n, s_inference[n], v_inference[n], s_inference[n] * v_inference[n]))
        final_result = sorted(final, reverse=True, key=lambda x: x[3])
        # print "final:\n",final_result
        with open(similarity_result_path+"result_" + str(inputpoints[0][3]) + "_" +  str(inputpoints[0][4])  + ".txt", 'wt') as fp:
            for i in final_result:
                fp.writelines(
                    str(i[0][0]) + ', ' + str(i[0][1]) + ',' + str(i[1]) + ', ' + str(i[2]) + ', ' + str(i[3]) + '\n')

        # endtime = time.time()
        # print endtime - starttime
        # final = sorted(final, reverse=True, key=lambda x: x[2])
        # if len(final) % 100 == 0:
        #     with open(similarity_result_path + brain, "wt")as ff:
        #         for neuron in final:
        #             ff.writelines(str(neuron[0]) + " " + str(neuron[1]) + " " + str(neuron[2]) + "\n")
        # with open(similarity_result_path + brain, "wt")as ff:
        #     for neuron in final:
        #         ff.writelines(str(neuron[0]) + " " + str(neuron[1]) + " " + str(neuron[2]) + "\n")
#
# thread_number=30
trunk_number=len(trunks)
thread1=threading.Thread(similarity_calculation(brain_dictionary,trunks[?*trunk_number/10:@*trunk_number/10],detect_list,cube_size))
# thread2=threading.Thread(similarity_calculation(brain_dictionary,trunks[1*trunk_number/thread_number:2*trunk_number/thread_number],detect_list,cube_size))
# thread3=threading.Thread(similarity_calculation(brain_dictionary,trunks[2*trunk_number/thread_number:3*trunk_number/thread_number],detect_list,cube_size))
# thread4=threading.Thread(similarity_calculation(brain_dictionary,trunks[3*trunk_number/thread_number:4*trunk_number/thread_number],detect_list,cube_size))
# thread5=threading.Thread(similarity_calculation(brain_dictionary,trunks[4*trunk_number/thread_number:5*trunk_number/thread_number],detect_list,cube_size))
# thread6=threading.Thread(similarity_calculation(brain_dictionary,trunks[5*trunk_number/thread_number:6*trunk_number/thread_number],detect_list,cube_size))
# thread7=threading.Thread(similarity_calculation(brain_dictionary,trunks[6*trunk_number/thread_number:7*trunk_number/thread_number],detect_list,cube_size))
# thread8=threading.Thread(similarity_calculation(brain_dictionary,trunks[7*trunk_number/thread_number:8*trunk_number/thread_number],detect_list,cube_size))
# thread9=threading.Thread(similarity_calculation(brain_dictionary,trunks[8*trunk_number/thread_number:9*trunk_number/thread_number],detect_list,cube_size))
# thread10=threading.Thread(similarity_calculation(brain_dictionary,trunks[9*trunk_number/thread_number:10*trunk_number/thread_number],detect_list,cube_size))
# thread11=threading.Thread(similarity_calculation(brain_dictionary,trunks[10*trunk_number/thread_number:11*trunk_number/thread_number],detect_list,cube_size))
# thread12=threading.Thread(similarity_calculation(brain_dictionary,trunks[11*trunk_number/thread_number:12*trunk_number/thread_number],detect_list,cube_size))
# thread13=threading.Thread(similarity_calculation(brain_dictionary,trunks[12*trunk_number/thread_number:13*trunk_number/thread_number],detect_list,cube_size))
# thread14=threading.Thread(similarity_calculation(brain_dictionary,trunks[13*trunk_number/thread_number:14*trunk_number/thread_number],detect_list,cube_size))
# thread15=threading.Thread(similarity_calculation(brain_dictionary,trunks[14*trunk_number/thread_number:15*trunk_number/thread_number],detect_list,cube_size))
# thread16=threading.Thread(similarity_calculation(brain_dictionary,trunks[15*trunk_number/thread_number:16*trunk_number/thread_number],detect_list,cube_size))
# thread17=threading.Thread(similarity_calculation(brain_dictionary,trunks[16*trunk_number/thread_number:17*trunk_number/thread_number],detect_list,cube_size))
# thread18=threading.Thread(similarity_calculation(brain_dictionary,trunks[17*trunk_number/thread_number:18*trunk_number/thread_number],detect_list,cube_size))
# thread19=threading.Thread(similarity_calculation(brain_dictionary,trunks[18*trunk_number/thread_number:19*trunk_number/thread_number],detect_list,cube_size))
# thread20=threading.Thread(similarity_calculation(brain_dictionary,trunks[19*trunk_number/thread_number:20*trunk_number/thread_number],detect_list,cube_size))
# thread21=threading.Thread(similarity_calculation(brain_dictionary,trunks[20*trunk_number/thread_number:21*trunk_number/thread_number],detect_list,cube_size))
# thread22=threading.Thread(similarity_calculation(brain_dictionary,trunks[21*trunk_number/thread_number:22*trunk_number/thread_number],detect_list,cube_size))
# thread23=threading.Thread(similarity_calculation(brain_dictionary,trunks[22*trunk_number/thread_number:23*trunk_number/thread_number],detect_list,cube_size))
# thread24=threading.Thread(similarity_calculation(brain_dictionary,trunks[23*trunk_number/thread_number:24*trunk_number/thread_number],detect_list,cube_size))
# thread25=threading.Thread(similarity_calculation(brain_dictionary,trunks[24*trunk_number/thread_number:25*trunk_number/thread_number],detect_list,cube_size))
# thread26=threading.Thread(similarity_calculation(brain_dictionary,trunks[25*trunk_number/thread_number:26*trunk_number/thread_number],detect_list,cube_size))
# thread27=threading.Thread(similarity_calculation(brain_dictionary,trunks[26*trunk_number/thread_number:27*trunk_number/thread_number],detect_list,cube_size))
# thread28=threading.Thread(similarity_calculation(brain_dictionary,trunks[27*trunk_number/thread_number:28*trunk_number/thread_number],detect_list,cube_size))
# thread29=threading.Thread(similarity_calculation(brain_dictionary,trunks[28*trunk_number/thread_number:29*trunk_number/thread_number],detect_list,cube_size))
# thread30=threading.Thread(similarity_calculation(brain_dictionary,trunks[29*trunk_number/thread_number:30*trunk_number/thread_number],detect_list,cube_size))
# similarity_calculation(brain_dictionary,trunks[0*trunk_number/4:0*trunk_number/4],detect_list,cube_size)
# Start
thread1.start()
# thread2.start()
# thread3.start()
# thread4.start()
# thread5.start()
# thread6.start()
# thread7.start()
# thread8.start()
# thread9.start()
# thread10.start()
# thread11.start()
# thread12.start()
# thread13.start()
# thread14.start()
# thread15.start()
# thread16.start()
# thread17.start()
# thread18.start()
# thread19.start()
# thread20.start()
# thread21.start()
# thread22.start()
# thread23.start()
# thread24.start()
# thread25.start()
# thread26.start()
# thread27.start()
# thread28.start()
# thread29.start()
# thread30.start()
# #
# #
# #
#
#