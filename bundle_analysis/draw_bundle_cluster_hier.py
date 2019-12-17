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


countpoints = 0
allpoints = OrderedDict()
exist = set()  # ***check if exist is needed to be outside the for loop
xyz=[]
time = 0
swc = " "
count=0
# ############somata
# x1 = []
# y1 = []
# z1 = []
# regex = re.compile(r"([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)")
# with open('soma_location_original.txt') as f:
#     for line in f:
#             # print line
#         for iter2 in regex.finditer(line):
#             x1.append(float(iter2.group(1))*10)
#             y1.append(float(iter2.group(2))*10)
#             z1.append(float(iter2.group(3))*10)
#
# fig = plt.figure()
# ax = fig.gca(projection='3d')
#
# plt.plot(x1, y1, z1, 'bo', Markersize=2)
#####################
# plt.show()
# plt.savefig("SO/"+i,dpi=300)
re_bundle=re.compile(
    r"([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\d]+)[\s]+([\d]+)"
)
# url=["cluster_10up.txt","cluster_100up.txt","cluster_50up.txt","cluster_200up.txt"]
# url=["static_cluster_10up.txt","static_cluster_3up.txt"]
url=["static_cluster_10up.txt"]

cluster_list=[]
neuron_collection=[]
for i in url:
    with open ("bundle_cluster_0p4/"+i,"rt") as fp:
        for j in fp:
            if j[-1]=='\n':
                cluster_list.append(j[0:-1])
            else:
                cluster_list.append(j)
for i in cluster_list:
    neuron_list=[]
    with open ("bundle_cluster_0p4/"+"sorting_"+i+".txt","rt") as fp:
        for j in fp:
            # print j
            if j[-1]=='\n':
                j=j[:-1]
                # k=j.split(",")
                name=j.split("_")
                neuron_list.append(str(nn.namelist.index(name[0]))+"_"+name[1])
            else:
                # k=j.split(",")
                name=j.split("_")
                neuron_list.append(str(nn.namelist.index(name[0]))+"_"+name[1])
    # print neuron_list
    neuron_collection.append(neuron_list)
R=[]
G=[]
B=[]
alpha=0.7
X=[]
Y=[]
Z=[]
SOMA_X=[]
SOMA_Y=[]
SOMA_Z=[]
color_index=1
fig = plt.figure()
ax = fig.gca(projection='3d')

for i in neuron_collection:
    color_change_stepsize = 0.03
    color_change = 0.1
    color_index_change=1/color_change_stepsize
    # soma_x=[]
    # soma_y=[]
    # soma_z=[]
    j = 0
    for neuron in i:
        x = []
        y = []
        z = []
        try:
            with open('bundles_c/' + neuron + ".txt") as f:  ##############原始資料的目錄
                j = 0
                for line in f:  # print line
                    for iter2 in re_bundle.finditer(line):
                        x.append(float(iter2.group(1)))
                        y.append(float(iter2.group(2)))
                        z.append(float(iter2.group(3)))
                        # if j == 0:
                        #     soma_x.append(float(iter2.group(1)))
                        #     soma_y.append(float(iter2.group(2)))
                        #     soma_z.append(float(iter2.group(3)))
                        #     j = j + 1
            if (color_index%(3*color_index_change)) > color_index_change*2:
                r = 1
                g = 1
                b = color_change_stepsize * (color_index%(3*color_index_change) - 2*color_index_change)
            elif (color_index%(3*color_index_change)) > color_index_change*1:
                r = 1
                g = color_change_stepsize * (color_index%(3*color_index_change) - 1*color_index_change)
                b = 0
            else:
                r = color_change_stepsize * (color_index%(3*color_index_change) - 0*color_index_change)
                g = 0
                b = 0
            color_index = color_index + 1
            plt.plot(x, y, z, "o", markersize=3, color=(r, g, b, 0.5))


            X.append(x)
            Y.append(y)
            Z.append(z)
        # SOMA_X.append(soma_x)
        # SOMA_Y.append(soma_y)
        # SOMA_Z.append(soma_z)
            R.append(r)
            G.append(g)
            B.append(b)
        except:
            print neuron
            continue
    for iii in xrange(10):
        for jjj in xrange(10):
            ax.view_init(iii * 36, jjj * 36)
            plt.savefig("bundle_cluster_0p4_fig/b_"+cluster_list[neuron_collection.index(i)]+"_"+str(iii)+"_"+str(jjj)+".png",dpi=100)

    # plt.show()
    plt.close()
    fig = plt.figure()
    ax = fig.gca(projection='3d')

for i in xrange(len(X)):
    print R[i],G[i],B[i]
    # plt.plot(SOMA_X[i], SOMA_Y[i], SOMA_Z[i], "o", markersize=10,color=(R[i]/10,G[i],B[i]/10,alpha))
    plt.plot(X[i], Y[i], Z[i], "o", markersize=3,color=(R[i],G[i],B[i],alpha))

for iii in xrange(10):
    for jjj in xrange(10):
        ax.view_init(iii * 36, jjj * 36)
        plt.savefig("bundle_cluster_0p4_fig/ALL_" + "_" + str(iii) + "_" + str(jjj) + ".png", dpi=100)
plt.close()
# plt.show()

#
#
# with open('cluster_1077.txt') as f: ##############要畫的cluster
#     for line in f:
#         # print line[:-1]
#         if line[len(line)-1]=='\n':
#             url.append(line[:-1])
#         else:
#             url.append(line)
# url2=[]
# with open('cluster_674.txt') as f:  ##############要畫的cluster
#     for line in f:
#         # print line[:-1]
#         if line[len(line) - 1] == '\n':
#             url2.append(line[:-1])
#         else:
#             url2.append(line)
# url = ["G0239-F-000001","G0239-F-000003","G0239-F-000004","G0239-F-000005","G0239-F-000006","G0239-F-000007","G0239-F-000008","G0239-F-000009" , "G0239-F-000010","G0239-F-000011","G0239-F-000012","G0239-F-000013"]
# SWC
# Insert form
# regex = re.compile(
#         r"([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\+|\-]?[\d]+)[\s]+([\d]+)"
#     )
# x2=[]
# y2=[]
# z2=[]
# soma_x2=[]
# soma_y2=[]
# soma_z2=[]
#
# for neuron in url2:
#     with open('Database_no modification/swc_spin/' + neuron+".swc") as f: ##############原始資料的目錄
#         j=0
#         for line in f:    # print line
#             for iter2 in regex.finditer(line):
#                 x2.append(float(iter2.group(3)))
#                 y2.append(float(iter2.group(4)))
#                 z2.append(float(iter2.group(5)))
#                 if j==0:
#                     soma_x2.append(float(iter2.group(3)))
#                     soma_y2.append(float(iter2.group(4)))
#                     soma_z2.append(float(iter2.group(5)))
#                     j=j+1

                # print iter2.group(2)


# # %%%%%%%%%%%%%%%%這邊開始是畫腦區%%%%%%%%%%%%%%%%%%%%%%%%%%
# # import new_merge as merge
# # fig = plt.figure()
# # ax = fig.gca(projection='3d')
# import pylab
#
#
# xyz = {}
#
#
# def radius(x, y, z, xm, ym, zm): #distance between two points
#     r = np.sqrt((x-xm)**2+(y-ym)**2+(z-zm)**2)
#     return r
#
# def prob(r, sigma): #probability formula
#     p = np.exp(-(r/sigma)**2/2)/(np.sqrt(2*(np.pi)**3)*(sigma)**3)
#     return p
#
# """
# def dotproduct(v1,v2):
#     return sum((a*b) for a, b in zip(v1, v2))
#
# def vecLength(v):
#     return np.sqrt(dotproduct(v,v))
#
# def angle(v1,v2):
#     return np.arccos(dotproduct(v1,v2))/ (vecLength(v1))*(vecLength(v2))
#
# def vecTheta(x1,y1,z1,x2,y2,z2):
#     theta = np.arccos((x1*x2+y1*y2+z1*z2)/radius(x=x1,y=y1,z=z1,xm=0,ym=0,zm=0)*radius(x=x2,y=y2,z=z2,xm=0,ym=0,zm=0))
#     return theta
# """ #previous code to compute vector angle
#
# def pos_to_int(x): #convert to int
#     return int(round(10*x))
#
# def unit_vec(vector): #convert to unit vector
#     return vector/ np.linalg.norm(vector)
#
# def angle_between(v1,v2): #compute angle betweeen vectors; use np.clip() to determine range; return in degrees
#     v1_u = unit_vec(v1)
#     v2_u = unit_vec(v2)
#     return 180*(np.arccos(np.clip(np.dot(v1_u,v2_u),-1.0,1.0)))/ np.pi
#
# def load_mesh(newurl,alpha): #function to draw neuronpils
#     allvertices = []
#     allfaces = []
#
#     for j in xrange(len(newurl)):
#         mesh = " " #can use "" as well
#         with open("Database_no modification/neuropils_mesh_txt/"+newurl[j]) as f: #open file that should be written into mesh
#             for line in f:
#                 mesh += line
#         vertices = []  # points
#         regex = re.compile(
#             r"^([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)[\s]+([\+|\-]?[\d]+[\.][\d]+)", re.M)
#         for iter in regex.finditer(mesh):
#             vertices.append([float(iter.group(1)), float(iter.group(2)), float(iter.group(3))])
#         allvertices.append(vertices)
#
#         faces = []  #constructed by every three points
#         regex = re.compile(
#             r"^([\d]+)[\s]+([\d]+)[\s]+([\d]+)", re.M)
#         for iter in regex.finditer(mesh):
#             faces.append([int(iter.group(1)), int(iter.group(2)), int(iter.group(3))])
#         allfaces.append(faces)
#
#     for j in xrange(len(allfaces)):
#         for i in allfaces[j]:
#             plt.plot([allvertices[j][i[0] - 1][0], allvertices[j][i[1] - 1][0]], [allvertices[j][i[0] - 1][1], allvertices[j][i[1] - 1][1]],
#                     [allvertices[j][i[0] - 1][2], allvertices[j][i[1] - 1][2]], 'r-', alpha=alpha)
#             plt.plot([allvertices[j][i[1] - 1][0], allvertices[j][i[2] - 1][0]], [allvertices[j][i[1] - 1][1], allvertices[j][i[2] - 1][1]],
#                     [allvertices[j][i[1] - 1][2], allvertices[j][i[2] - 1][2]], 'r-', alpha=alpha)
#             plt.plot([allvertices[j][i[2] - 1][0], allvertices[j][i[0] - 1][0]], [allvertices[j][i[2] - 1][1], allvertices[j][i[0] - 1][1]],
#                     [allvertices[j][i[2] - 1][2], allvertices[j][i[0] - 1][2]], 'r-', alpha=alpha)



# newurl = ["eb_8_instd_l.txt"] # ,"lh_25_instd_l.txt","mb_4_instd_l.txt"
# load_mesh(newurl=newurl,alpha=0.05)
# plt.rcParams['savefig.dpi'] = 500 #set the resolution for fig




#################################################################

#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.add_collection3d

#for key,value in x.allpoints.iteritems():
#    print  key[0], key[1], key[2],x.allpoints.get(key)[0],value #, value[1], value[2], value[3]key[0], key[1], key[2]
    #if value != None:
    #    ax.quiver(key[0], key[1], key[2], value[1], value[2], value[3], length=0.1)


#print allpoints

# for i in merge.mergedCube:
#     ax.quiver(i[0], i[1], i[2], merge.mergedCube[i][0], merge.mergedCube[i][1], merge.mergedCube[i][2], length=200,alpha=.2,arrow_length_ratio=.3,linewidth=1.5)


# for ii in xrange(0,180,20): #draw the fig every 20 degrees in azumuth
#     ax.view_init(elev=0., azim=ii)
#     pylab.savefig("al&mb_al&lh_azimuth%d_w15.png" %ii)

# ax.view_init(elev=90, azim=0) #draw the fig in top view
# pylab.savefig("al&mb_al&lh_top_view_w15.png" , bbox_inches='tight')


# ax.view_init(azim=0, elev=90) #set the view point
# plt.show()




