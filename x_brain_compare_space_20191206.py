#! /usr/bin/env python
# coding: utf-8
import numpy as np
import re
import time
from collections import OrderedDict
import cPickle as pickle
import x_infer_spatial_space_20191206 as s_infer
import x_infer_vec_space_20191206 as v_infer
import os
import random
import bound_list_sig20um_2sig as bound  # sphere bound
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

p_threshold="19000/"
brain_pathway=p_threshold+"19000_swc_warped_no_replicate_fragmentize_brains/"
space_pathway=p_threshold+"19000_swc_warped_no_replicate_fragmentize_space_100_pickle/"
result_pathway=p_threshold+"19000_swc_warped_no_replicate_fragmentize_bundle_result_20/"
tmp_pathway=p_threshold+"19000_swc_warped_no_replicate_fragmentize_bundle_result_tmp_20/"

if not os.path.isdir(tmp_pathway):
    os.mkdir(tmp_pathway)
if not os.path.isdir(result_pathway):
    os.mkdir(result_pathway)


with open ("19000_brain_list.txt","rt") as f:
    for line in f:
        url.append(line[:-1])
print len(url)
print url

url=os.paht.listdir(brain_pathway)

c=0
regex = re.compile(
    # ptr_x,ptr_y,ptr_z,name_index,branch_id,total_node,0,0,0,son_number
    r"([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\d]+)"
)

#10600_77.txt---短版(20)的暫存
# url_index=url.index("103_206.txt")
# url[url_index:]

###Build detect
detect_list = []  # list: the sphere bound for the input point→建一個以input point為中心的球狀座標表單
for j in bound.bounding_list:
    detect_list.append([pos_to_int(j['x']),pos_to_int(j['y']),pos_to_int(j['z']),math.log10(j['p'] / minp) / stdp])
detect_list = sorted(detect_list, reverse=True, key=lambda x: x[3])
print detect_list[0]


length_list=len(url)
url_index=0

##Check to avoid repeat
#Dictionary=[]
#for dirPath, dirNames, fileNames in os.walk(result_pathway):
#    for i in fileNames:
#        grouping=i.split("_")
#        if int(grouping[1]) not in Dictionary:
#            Dictionary.append(int(grouping[1]))

for name in url:
    print name
    #Here tmp_list is for brain~~~I will check branch id by finish the branch
    tmp_list=os.path.listdir(tmp_pathway)
    branch_list=[]
    if name in tmp_list:
        print "exist:", name
        branch_list=os.path.listdir(result_pathway+name[:-4])
        continue  ##if check branch~delete this line
        
    else:
        with open(tmp_pathway+name,"wt")as ff:
            ff.write("exist")
        os.mkdir(result_pathway+name[:-4])
    url_index=url_index+1
    length = 0
    # name=url[sampling_count]
    starttime=time.time()
    swc = " "
    with open(brain_pathway+name) as f:
        pre_branch_id = 0
        length=0
        inputpoints = []
        for line in f:
            for iter in regex.finditer(line):
                branch_id = int(iter.group(5))
                if "result_" + name[:-4] + "_" + str(pre_branch_id) + "_length_" + str(length) + ".txt" in branch_list:
                    print "exist:", name[:-4]+"_"+str(pre_branch_id)
                    continue
                else:
                    if branch_id != pre_branch_id:
                        if length < 20:
                            pre_branch_id = branch_id
                            length = 0
                            del inputpoints
                            inputpoints = []
                        else:
                            inference = s_infer.inference(inputpoints=inputpoints,detect_list=detect_list)
                            inference_vec = v_infer.inference_vec(inputpoints=inputpoints,detect_list=detect_list)
                            print "succeed"
                            print name+str(branch_id)
                            #print "inference=", inference
                            #print "inference_vec=", inference_vec
                            final = []
                            for n in inference:
                                if n not in inference_vec:
                                    inference_vec[n] = 0
                                final.append((n, inference[n], inference_vec[n], inference[n] * inference_vec[n]))
                            final_result = sorted(final, reverse=True, key=lambda x: x[3])
                            # print "final:\n",final_result
                            with open(result_pathway+name[:-4]+"/result_" + name[:-4] + "_" + str(pre_branch_id) + "_length_" + str(length) + ".txt", 'w') as fp:
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
                    else:
                        inputpoints.append([int(iter.group(1)), int(iter.group(2)), int(iter.group(3)), int(iter.group(7)),
                                        int(iter.group(8)), int(iter.group(9))])
                        length = length + 1
    if length >20 :
        inference = s_infer.inference(inputpoints=inputpoints,detect_list=detect_list)
        inference_vec = v_infer.inference_vec(inputpoints=inputpoints,detect_list=detect_list)
        print "inference=", inference
        print "inference_vec=", inference_vec
        final = []
        for n in inference:
            if n not in inference_vec:
                inference_vec[n] = 0
            final.append((n, inference[n], inference_vec[n], inference[n] * inference_vec[n]))
        final_result = sorted(final, reverse=True, key=lambda x: x[3])
        # print "final:\n",final_result
        with open(result_pathway+name[:-4]+"/result_" + name[:-4] + "_" + str(
                branch_id) + "_length_" + str(length) + ".txt", 'w') as fp:
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
