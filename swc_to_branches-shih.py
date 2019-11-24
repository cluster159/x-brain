#! /usr/bin/env python
# coding: utf-8
import requests
import re
import os
import numpy as np
import time

##Variable: threshold,p
##############################################################################
threshold="19000"
p="x_brain/"
pathway_original=p+"brain_clean/tracing_bin1/"+threshold
# pathway_transformed=p+threshold+"_branches_transformed"
pathway_unwarped=p+threshold+"_branches_unwarped/"
pathway_warped=p+threshold+"_branches_warped/"
pathway_final=p+threshold+"_branches_final/"
pathway_swc_warped=p+threshold+"_swc_warped/"
pathway_swc_unwarped=p+threshold+"_swc_unwarped/"
# ############################################################################

#######################################################
re_hx=re.compile(
    r"([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)"
)

re_branches=re.compile(
    r"([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)"
)

regex = re.compile(
    # ptr_x,ptr_y,ptr_z,name_index,branch_id,total_node,0,0,0,SO,son_number
    r"([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\d]+)[\s]+([\d]+)"
)
re_pair_swc=re.compile(
    r"([\d]+[\.][\d]+)[\s]+([\d]+[\.][\d]+)[\s]+([\d]+[\.][\d]+)"
)
##########################################################


#########get warping list############
warping_list = []
for dirPath, dirNames, fileNames in os.walk(pathway_swc_warped):
    for i in fileNames:
        warping_list.append(i[:-21])
print warping_list
print len(warping_list)

########先進行transform
try:
    # os.makedirs(pathway_transformed)
    os.makedirs(pathway_unwarped)
    os.makedirs(pathway_warped)
    os.makedirs(pathway_final)
except:
    print "error"

url=[]
list=[]
hx=[]
for dirPath, dirNames, fileNames in os.walk(pathway_original):
    for i in fileNames:
        if i.find("Branches.txt")!=-1:
            print i
            if i[:-19] in warping_list:
                list.append(i)
                url.append(dirPath.replace('\\', '/'))
        if i.find("Script.hx")!=-1:
            if i [:-16] in warping_list:
            # print i[:-16]
                hx.append(i)
print list
print hx
print len(list)


for i in xrange(len(list)):
    swc_warped_X = []
    swc_warped_Y = []
    swc_warped_Z = []

    swc_unwarped_X = []
    swc_unwarped_Y = []
    swc_unwarped_Z = []

    unwarped_size = 0
    warped_size = 0
    check_length = 0
    with open(url[i]+"/"+hx[i],"r") as f:
        x_transform = []
        y_transform = []
        z_transform = []
        R_transform = []
        for line in f:
            if line.find("setTransform")!=-1:
                for iter in re_hx.finditer(line):
                    x_transform.append(float(iter.group(1)))
                    x_transform.append(float(iter.group(5)))
                    x_transform.append(float(iter.group(9)))
                    x_transform.append(float(iter.group(13)))
                    y_transform.append(float(iter.group(2)))
                    y_transform.append(float(iter.group(6)))
                    y_transform.append(float(iter.group(10)))
                    y_transform.append(float(iter.group(14)))
                    z_transform.append(float(iter.group(3)))
                    z_transform.append(float(iter.group(7)))
                    z_transform.append(float(iter.group(11)))
                    z_transform.append(float(iter.group(15)))
                    R_transform.append(float(iter.group(4)))
                    R_transform.append(float(iter.group(8)))
                    R_transform.append(float(iter.group(12)))
                    R_transform.append(float(iter.group(16)))
                    # print "TM"
                    # print x_transform
                    # print y_transform
                    # print z_transform
                    # print R_transform
                break
    with open(url[i]+"/"+list[i],"r") as f:
        X=[]
        Y=[]
        Z=[]
        R=[]
        for line in f:
            for iter in re_branches.finditer(line):
                x=float(iter.group(3))
                y=float(iter.group(2))
                z=float(iter.group(1))
                r=1
                x_t=x*x_transform[0]+y*x_transform[1]+z*x_transform[2]+r*x_transform[3]
                y_t=x*y_transform[0]+y*y_transform[1]+z*y_transform[2]+r*y_transform[3]
                z_t=x*z_transform[0]+y*z_transform[1]+z*z_transform[2]+r*z_transform[3]
                r_t=x*R_transform[0]+y*R_transform[1]+z*R_transform[2]+r*R_transform[3]
                X.append(x_t)
                Y.append(y_t)
                Z.append(z_t)
                R.append(r_t)
                # print "original: ", x, y, z
                # print "transform: ", x_t,y_t,z_t
                check_length=check_length+1
    unwarped_size = len(X)
    path = pathway_unwarped + list[i][:-4] + "_before_warp.txt"
    check_write_length = 0
    with open(path,"w") as f:
        for j in xrange(len(X)):
            f.writelines(str(X[j])+" "+str(Y[j])+" "+str(Z[j])+"\n")
            check_write_length = check_write_length + 1
    print path
    print len(X)
    print "branches_extract_transform_over"
######################construct warping dictionary ####################
    with open(pathway_swc_warped+warping_list[i]+"-brain_after_warp.txt","rt")as f_swc:
        for line in f_swc:
            for iter in re_pair_swc.finditer(line):
                swc_warped_X.append(float(iter.group(1)))
                swc_warped_Y.append(float(iter.group(2)))
                swc_warped_Z.append(float(iter.group(3)))
    with open(pathway_swc_unwarped+warping_list[i]+"-brain_before_warp.txt","rt")as f_swc:
        for line in f_swc:
            for iter in re_pair_swc.finditer(line):
                swc_unwarped_X.append(float(iter.group(1)))
                swc_unwarped_Y.append(float(iter.group(2)))
                swc_unwarped_Z.append(float(iter.group(3)))
    print pathway_swc_unwarped+warping_list[i]+"-brain_before_warp.txt"
    print pathway_swc_warped+warping_list[i]+"-brain_after_warp.txt"
    print list[i]
    for pair_index in xrange(len(X)):
        for swc_index in xrange(len(swc_unwarped_X)):
            check_check=0
            if X[pair_index]==swc_unwarped_X[swc_index] and Y[pair_index]==swc_unwarped_Y[swc_index] and Z[pair_index]==swc_unwarped_Z[swc_index]:
                X[pair_index] = swc_warped_X[swc_index]
                Y[pair_index] = swc_warped_Y[swc_index]
                Z[pair_index] = swc_warped_Z[swc_index]
                check_check=1
                break
        print len(X), pair_index, pair_index/len(X)
        if check_check==0:
            print "fail"
            print pair_index
            print X[pair_index],Y[pair_index],Z[pair_index]
            break

    with open(pathway_warped + list[i] + "_after_warp.txt", 'wt')as fout:
        for pair_index in len(X):
            fout.writelines(str(X[pair_index])+" "+str(Y[pair_index])+" "+str(Z[pair_index])+"\n")
        print "warping_over"

################################################Check and return to original form###################
    f_name = list[i][:-4]
    with open(pathway_final + f_name + "_shih.txt", "w") as fp:
        ptr = 0
        with open(url[i] + "/" + list[i], "r") as fo:
            for line in fo:
                # print line
                check = 0
                for iter in re_branches.finditer(line):
                    check = 1
                    r1 = iter.group(4)
                    r2 = iter.group(5)
                    r3 = iter.group(6)
                    r4 = iter.group(7)
                if check == 1:
                    fp.writelines(str(Z[ptr]) + " " + str(Y[ptr]) + " " + str(
                        X[ptr]) + " " + r1 + " " + r2 + " " + r3 + " " + r4 + "\n")
                    ptr = ptr + 1
                else:
                    # print line
                    fp.writelines(line)
    print "final_confirmation"
    free(swc_unwarped)
    free(swc_warped)


