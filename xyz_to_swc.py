#! /usr/bin/env python
# coding: utf-8
import re
import os
import numpy as np
import time

##Variable: threshold,p
#########get warping list############

##########################################
regex = re.compile(
    # ptr_x,ptr_y,ptr_z,name_index,branch_id,total_node,0,0,0,SO,son_number
    r"([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\d]+)[\s]+([\d]+)"
)
re_warp=re.compile(
    r"([\+|\-]?[\d]+[\.]?[\d]?)[\s]+([\+|\-]?[\d]+[\.]?[\d]?)[\s]+([\+|\-]?[\d]+[\.]?[\d]?)"
)

re_swc=re.compile(
    r"([\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+)"
)


##########################################
threshold="19000"
p="E:/20191115_new_bin1_xbrain/"
pathway_original=p+threshold+"_swc_original/"
pathway_warped=p+threshold+"_xyz_warped/"
#C:\Fly brain database\x2fc_affine\v1\20000_swc_unwarped
pathway_final=p+threshold+"_swc_warped/"
try:
    # os.makedirs(pathway_unwarped)
    # os.makedirs(pathway_warped)
    os.makedirs(pathway_final)
except:
    print "exist"

url=[]
list=[]
for dirPath, dirNames, fileNames in os.walk(pathway_warped):
    for i in fileNames:
        if i.find("FC")!=-1:
            list.append(i)

#-38
for i in xrange(len(list)):
    ################################################Check and return to original form###################
    f_name=list[i][:-16]
    print f_name
    with open(pathway_original+list[i][:-16] + ".swc", "r") as fo:
        r1 = []
        r2 = []
        r3 = []
        r4 = []
        for line in fo:
            for iter in re_swc.finditer(line):
                r1.append(iter.group(1))
                r2.append(iter.group(2))
                r3.append(iter.group(6))
                r4.append(iter.group(7))
    original_size=len(r4)
    X_w=[]
    Y_w=[]
    Z_w=[]
    with open(pathway_warped + list[i], "rt") as fp:
        count=0
        start=0
        for line in fp:
            if line =="\n":
                break
            if line[0].find("@")!=-1:
                start=1
                continue
            if start==1:
                group = line.split(" ")
                if len(group) < 3:
                    print "warping error"
                    # break
                # print group
                X_w.append(float(group[0]))
                Y_w.append(float(group[1]))
                if group[2][-1] == "\n":
                    group[2] = group[2][:-1]
                Z_w.append(float(group[2]))
    warped_size = len(X_w)
    print "check: ",original_size, warped_size
    with open(pathway_final + f_name + "_x2fc.swc", "wt") as fp:
        for ptr in xrange(len(r1)):
            fp.writelines(r1[ptr] + " " + r2[ptr] + " " +str(X_w[ptr]) + " " + str(Y_w[ptr]) + " " + str(Z_w[ptr]) + " " + r3[ptr] + " " + r4[ptr] + "\n")
    print "final_confirmation"




