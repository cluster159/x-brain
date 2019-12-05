#! /usr/bin/env python
# coding: utf-8
import numpy as np
import re
import time
from collections import OrderedDict
import cPickle as pickle
import x_infer_spatial_space as s_infer
import x_infer_vec_space as v_infer
import os
import random
start = time.time()

class Scalar:
    value = 0.5
    sigma = 2.0
def radius(x, y, z, xm, ym, zm): #distance between two points
    r = np.sqrt((x-xm)**2+(y-ym)**2+(z-zm)**2)
    return r

def pos_to_int(x): #convert to int(x*10)
    return int(round(10*x))

def name_to_index(name): #the input name should be in "filename" form
    return nn.namelist.index(name)

def toScalar(x): #convert to resolution
    if round((x%Scalar.value)/Scalar.value) == 1:
        x = x-x%Scalar.value+Scalar.value
    else:
        x = x-x%Scalar.value
    return float("%.3f" %x)


names=[]
url=[]

with open ("x_swc_warped_20000_list.txt","rt") as f:
    for line in f:
        url.append(line[:-1])
print len(url)
print url

c=0
regex = re.compile(
    # ptr_x,ptr_y,ptr_z,name_index,branch_id,total_node,0,0,0,son_number
    r"([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\d]+)"
)

#10600_77.txt---短版(20)的暫存
# url_index=url.index("103_206.txt")
# url[url_index:]
length_list=len(url)
url_index=0
Dictionary=[]
for dirPath, dirNames, fileNames in os.walk("x_brain_result_long_20_20000"):
    for i in fileNames:
        grouping=i.split("_")
        if int(grouping[1]) not in Dictionary:
            Dictionary.append(int(grouping[1]))

for name in url:
    check_name=name
    grouping=check_name.split(".")
    if int(grouping[0]) in Dictionary:
        print "brain exist", int(grouping[0])
        continue
    else:
        print "current:", grouping[0]

    url_index=url_index+1
    length = 0
    # name=url[sampling_count]
    starttime=time.time()
    print name
    swc = " "
    with open("x_brain_c_0p3_20000/"+name) as f:
        pre_branch_id = 0
        length=0
        inputpoints = []
        for line in f:
            for iter in regex.finditer(line):
                branch_id = int(iter.group(5))
                if os.path.isfile("x_brain_result_long_20_20000_"+str(length_list/31)+"/result_" + name[:-4] + "_" + str(
                                    pre_branch_id) + "_length_" + str(length) + ".txt"):
                    print "exist:", name[:-4]+"_"+str(branch_id)
                    continue
                else:
                    if branch_id != pre_branch_id:
                        #print "Here", url_index, "x_brain_result_long_50_20000/result_" + name[:-4] + "_" + str(
                                   # branch_id) + "_length_" + str(length) + ".txt"
                        if length < 20:
                            pre_branch_id = branch_id
                            length = 0
                            del inputpoints
                            inputpoints = []
                        else:
                            inference = s_infer.inference(inputpoints=inputpoints)
                            inference_vec = v_infer.inference_vec(inputpoints=inputpoints)
                            print "succeed"
                            print "Here", url_index, "x_brain_result_long_20_20000_"+str(length_list/31)+"/result_" + name[:-4] + "_" + str(
                                    pre_branch_id) + "_length_" + str(length) + ".txt"

                            #print "inference=", inference
                            #print "inference_vec=", inference_vec
                            final = []
                            for n in inference:
                                if n not in inference_vec:
                                    inference_vec[n] = 0
                                final.append((n, inference[n], inference_vec[n], inference[n] * inference_vec[n]))
                            final_result = sorted(final, reverse=True, key=lambda x: x[3])
                            # print "final:\n",final_result
                            with open("x_brain_result_long_20_20000_"+str(length_list/31)+"/result_" + name[:-4] + "_" + str(
                                    pre_branch_id) + "_length_" + str(length) + ".txt", 'w') as fp:
                                for i in final_result:
                                    fp.writelines(
                                        str(i[0][0]) + ', ' + str(i[0][1]) + ',' + str(i[1]) + ', ' + str(
                                            i[2]) + ', ' + str(i[3]) + '\n')
                            endtime = time.time()
                            print endtime - starttime
                            pre_branch_id = branch_id
                            length = 0
                            del inputpoints
                            inputpoints = []
                    inputpoints.append([int(iter.group(1)), int(iter.group(2)), int(iter.group(3)), int(iter.group(7)),
                                        int(iter.group(8)), int(iter.group(9))])
                    length = length + 1
    if length >20 :
        inference = s_infer.inference(inputpoints=inputpoints)
        inference_vec = v_infer.inference_vec(inputpoints=inputpoints)
        print "inference=", inference
        print "inference_vec=", inference_vec
        final = []
        for n in inference:
            if n not in inference_vec:
                inference_vec[n] = 0
            final.append((n, inference[n], inference_vec[n], inference[n] * inference_vec[n]))
        final_result = sorted(final, reverse=True, key=lambda x: x[3])
        # print "final:\n",final_result
        with open("x_brain_result_long_20_20000_"+str(length_list/31)+"/result_" + name[:-4] + "_" + str(
                pre_branch_id) + "_length_" + str(length) + ".txt", 'w') as fp:
            for i in final_result:
                fp.writelines(
                    str(i[0][0]) + ', ' + str(i[0][1]) + ',' + str(i[1]) + ', ' + str(i[2]) + ', ' + str(
                        i[3]) + '\n')
        endtime = time.time()
        print endtime - starttime
        pre_branch_id = branch_id
        length = 0
        del inputpoints
        inputpoints = []
