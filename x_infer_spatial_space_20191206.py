#! /usr/bin/env python
# coding: utf-8

#如何處理branch??????
##要改的部分:
# 1. inference直接改成相乘，每次做一次normalize
# 2. 要內插→直接從原本的xyz去輸出即可!!!!!
# 3. 只考慮最高的
# 4. 新加入的權重給漢原本排名中最低的一樣(可調整)
# 5. 要做100個data給老師
# 6. 要考慮branch
# 7. 消失的權重給0
# 8. 重新出現的話?→先跟4一樣

#暫時不考慮單純的vector(改寫成兩者相乘的)

# 已確認內切幾乎一樣(除了浮點誤差)
# 已確認分割的點沒問題
# 已確認sigma不同大小不影響最終的分數(500的分數應該是對的)(之後在確認一次)但和125的不同!!!WHY?


import numpy as np
import math as math
#import pandas as pd
import re
#from itertools import cycle
#import operator
import time
#import sphere as s
from collections import OrderedDict
# from matplotlib import pyplot as plt
import cPickle as pickle
# fig = plt.figure()
cube_size=100
cube_number=20
p_threshold="19000/"
brain_pathway=p_threshold+"19000_swc_warped_no_replicate_fragmentize_brains/"
space_pathway=p_threshold+"19000_swc_warped_no_replicate_fragmentize_space_100_pickle/"
result_pathway=p_threshold+"19000_swc_warped_no_replicate_fragmentize_bundle_result_50/"
tmp_pathway=p_threshold+"19000_swc_warped_no_replicate_fragmentize_bundle_result_tmp_50/"

xyz = {}
class Scalar: #some fixed parameters
    value = 0.5 #value of resolution
    sigma = 2.0 #sigma that can effect the size of sphere bound

# 注意!!!minp要隨著bounding list變動，另外小數位數要取對

neuro_space=[]
neuro_space_list=[]

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
    return float("%.3f" %x)

def name_to_index(name): #convert name to int
    return nn.namelist.index(name)

def pos_to_int(x): #convert to int(x*10)
    return int(round(10*x))

def theta_to_int(x): #convert to int(x*100)
    return int(round(100*x))

def unit_vec(vector):
    return vector/ np.linalg.norm(vector)

def angle_between(v1,v2): #use np.clip() to avoid error; return value in radians
    v1_u = unit_vec(v1)
    v2_u = unit_vec(v2)
    return np.arccos(np.clip(np.dot(v1_u,v2_u),-1.0,1.0))

def neurite_vec(startpoint,endpoint): ##??
    return (endpoint[0]-startpoint[0],endpoint[1]-startpoint[1],endpoint[2]-startpoint[2])
#endpoint-startpoint

def add_neuro_space(mark_x, xindex, mark_y, yindex,mark_z, zindex):
    # print "add"
    neuro_space_name = "space_x"
    if mark_x >= 0:
        neuro_space_name += "p"
        neuro_space_name += str(xindex)
    elif mark_x < 0:
        neuro_space_name += "n"
        neuro_space_name += str(xindex)

    neuro_space_name += "y"
    if mark_y >= 0:
        neuro_space_name += "p"
        neuro_space_name += str(yindex)
    elif mark_y < 0:
        neuro_space_name += "n"
        neuro_space_name += str(yindex)

    neuro_space_name += "z"
    if mark_z >= 0:
        neuro_space_name += "p"
        neuro_space_name += str(zindex)
    elif mark_z < 0:
        neuro_space_name += "n"
        neuro_space_name += str(zindex)
    neuro_space_name += ".txt.pickle"
    try:
        with open(space_pathway+neuro_space_name) as fp:
            datapoints = pickle.load(fp)
            #print "QQQQQQ",datapoints
            neuro_space.append(datapoints)
            neuro_space_list.append((mark_x,xindex,mark_y, yindex,mark_z, zindex))
        fp.close()
    except:
#        with open("xbrain_open_error_100.txt","a") as fp:
 #           fp.writelines('('+str(xindex)+','+str(yindex)+','+str(zindex)+')')
  #      fp.close()
        pass

def detect_neuron(inputx, inputy, inputz, detect_list): #version for input is a single point; datapoints is an OrderedDict
    t1=time.time()
    #print "HAHAHA",inputx, inputy, inputz
    involved_neurons = OrderedDict() #record neuron with its possibility ((key:n名稱):(value: p機率),(n):(p),...)
    for i in detect_list:
        list_index = -1
        x=i['x']+inputx
        y=i['y']+inputy
        z=i['z']+inputz
        xindex = abs(x) / cube_size
        if x >= 0:
            mark_x=1
        elif x < 0:
            mark_x=-1
        yindex = abs(y) / cube_size
        if y >= 0:
            mark_y=1
        elif y < 0:
            mark_y=-1
        zindex = abs(z) / cube_size
        if z >= 0:
            mark_z=1
        elif z < 0:
            mark_z=-1
        if (mark_x, xindex, mark_y, yindex,mark_z, zindex) not in neuro_space_list:
            add_neuro_space(mark_x, xindex, mark_y, yindex,mark_z, zindex)
            for kk in neuro_space_list:
                if abs(kk[1]*kk[0] - xindex*mark_x) > cube_number or abs(kk[2]*kk[3] - mark_y*yindex) > cube_number or abs(kk[4]*kk[5] - mark_z*zindex) > cube_number:
                    # print neuro_space
                    neuro_space.remove(neuro_space[neuro_space_list.index(kk)])
                    neuro_space_list.remove(kk)
        try:
            list_index=neuro_space_list.index((mark_x,xindex,mark_y,yindex,mark_z,zindex))
        except:
            pass
        if list_index>=0:
            if (x,y,z) in neuro_space[list_index]:
                for k in neuro_space[list_index][(x,y,z)]:
                    if (k[0],k[1]) not in involved_neurons:
                        involved_neurons[(k[0],k[1])] = i['p']
                    else:
                        if involved_neurons[(k[0],k[1])] < i['p']:
                            involved_neurons[(k[0],k[1])] = i['p']
    t6=time.time()
    # print "total", t6-t1
    return involved_neurons #OrderedDict



def inference(inputpoints,detect_list): #assume inputpoints is a list of lists/tuples; datapoints is an OrderedDict，這個沒包含vector
    #print "infer"
    neuron_number = 0
    global result
    result = OrderedDict()
    count = 0
    for i in xrange(len(inputpoints)):
        # print "ptr",i
        involved_points = detect_neuron(inputx=inputpoints[i][0],inputy=inputpoints[i][1],inputz=inputpoints[i][2], detect_list=detect_list)
        # print 'check',involved_points
        neuron_new_number=len(involved_points)
        #involved_neurons is an OrderedDict((n):(p),(n):(p)...)
        if neuron_new_number == 0 and neuron_number ==0:
            continue
        else:
            if neuron_new_number>neuron_number:
                neuron_number = neuron_new_number
            count = count + 1
            if count == 1:
                result = involved_points
            else:
                for m in involved_points:
                    if m in result:
                        # print 'm=',m
                        result[m]=result[m]+involved_points[m]
                    else:
                        result[m] = involved_points[m]

    #             四種狀況，在result(舊的)和involved(新的)都出現，只出現在result，只出現在involved，兩個都沒出現過
    #             都出現的不管是哪個人跑迴圈都可以偵測到，只出現在result的，不需做處理，因為直接加上0，而只出現在inovlved，則需要加上新的值(過去的都自動設為0，所以都不用處理)，如果是兩個都沒出現，當然不用處理
    #              →因此只需要跑involved的迴圈就可!!!!!!!!!!
    # print 'number',neuron_new_number,neuron_number
    # print 'all',result
    # del neuro_space_list
    # del neuro_space
    neuron_length=len(inputpoints)
    for i in result:
        result[i]=result[i]/neuron_length
    return result


def inference_step_by_step(inputpoints,url): #assume inputpoints is a list of lists/tuples; datapoints is an OrderedDict
    # url = ['Cha-F-100237', 'TH-F-000099', 'TH-F-000091', 'TH-F-100083', 'TH-F-100095']
    neuron_number=0
    node_count=0
    neuron_index_startpoint=OrderedDict()
    detect_list = []  # list: the sphere bound for the input point
    for j in bound.bounding_list:
        xyz = {
            'x': pos_to_int(j['x']),
            'y': pos_to_int(j['y']),
            'z': pos_to_int(j['z']),
            'p': math.log10(j['p'] / minp)/math.log10(maxp/minp)
        }
        xyz = xyz.copy()
        detect_list.append(xyz)
    global result
    result = OrderedDict()
    result_count = OrderedDict()
    result_average = OrderedDict()
    result_step_by_step = []
    result_step_by_step_count = []
    result_step_by_step_average = []
    for i in xrange(0,len(inputpoints)):
        node_count=node_count+1
        involved_points = detect_neuron(inputx=inputpoints[i][0], inputy=inputpoints[i][1],inputz=inputpoints[i][2], detect_list=detect_list)
        neuron_new_number = len(involved_points)
        # print "10287:",involved_points[10287]
        # print len(involved_points),involved_points
        # print involved_points
        # involved_neurons is an OrderedDict((n):(p),(n):(p)...)
        if neuron_number == 0 and neuron_new_number==0:
            continue
        else:
            if neuron_number<neuron_new_number:
                neuron_number=neuron_new_number
            for m in involved_points:
                if m in neuron_index_startpoint:
                        # print 'm=',m
                    result[m]=result[m]+involved_points[m]
                    result_count[m]=result_count[m]+1
                    # print 'In'
                else:
                    neuron_index_startpoint[m]=i
                    # print neuron_index_startpoint
                    result[m] = involved_points[m]
                    result_count[m]=1
            for n in result:
                tmp =result[n]/node_count
                result_average[n]=tmp
                # print n
                # print result[n], result_count[n], result_average[n]
            result_step_by_step.append(result)
            result_step_by_step_count.append(result_count)
            result_step_by_step_average.append(result_average)
            result=result.copy()
            result_count=result_count.copy()
            result_average=result_average.copy()
    # print result_step_by_step_average
    # print 'neuron_number_original',neuron_number
    node_number=0
    neuron_number = 0
    neuron_index=[]
    node_theta=[]
    node=[]
    neuron=[]
    neuron_theta=[]
    node_everynode=[]
    node_theta_everynode=[]
    neuron_everynode=[]
    neuron_theta_everynode=[]
    for m in neuron_index_startpoint:
        neuron_index.append(m)
        for n in xrange(neuron_index_startpoint[m],len(result_step_by_step)):
            node_theta_everynode.append(result_step_by_step[n][m])
            node_everynode.append(n)
            node_theta.append(result_step_by_step_average[n][m])
            node.append(n)
            node_number = node_number + 1
        neuron_everynode.append(node_everynode)
        neuron_theta_everynode.append(node_theta_everynode)
        neuron.append(node)
        neuron_theta.append(node_theta)
        node = []
        node_theta = []
        node_everynode = []
        node_theta_everynode = []
        neuron_number = neuron_number + 1
    # print "here"
    # print neuron_everynode
    # print neuron_theta_everynode
    print "neuron_number",len(neuron)
    # print neuron

    fig, ax = plt.subplots()
    for m in xrange(len(neuron)):
        # if nn.namelist[neuron_index[m]] not in url:
        #     continue
        ax.plot(neuron[m], neuron_theta[m],label=str(nn.namelist[neuron_index[m]]))
        ax.legend(loc='upper right',fontsize='xx-large')
        ax.set_ylim([0, 1.1])
        ax.set_xlim([0, len(result_step_by_step)+100])
        ax.set_title("Mean Score",fontsize=40)
        ax.set_ylabel("Score",fontsize=40)
        ax.set_xlabel("Node number",fontsize=40)
        ax.tick_params(labelsize=33)
    plt.show()

    fig, ax = plt.subplots()
    for m in xrange(len(neuron)):
        # fig, ax = plt.subplots()
        # if nn.namelist[neuron_index[m]] not in url:
        #     continue
        ax.plot(neuron_everynode[m], neuron_theta_everynode[m],label=str(nn.namelist[neuron_index[m]]))
        ax.legend(loc='upper left',fontsize='xx-large')
        # ax.set_ylim([0, 2500])
        ax.set_xlim([0, len(result_step_by_step)+100])
        ax.set_title("Score", fontsize=40)
        ax.set_ylabel("Score", fontsize=40)
        ax.set_xlabel("Node number", fontsize=40)
        ax.tick_params(labelsize=33)
    plt.show()

    for m in xrange(len(neuron)):
        print m
        print neuron_index[m]
        print neuron[m]
        print neuron_theta[m]
        fig, ax = plt.subplots()
        ax.plot(neuron[m],neuron_theta[m],label=str(nn.namelist[neuron_index[m]]))
        ax.set_ylim([0, 1.1])
        ax.set_xlim([0, len(result_step_by_step)+100])
        ax.set_title("Mean score", fontsize=25)
        ax.set_ylabel("Score", fontsize=25)
        ax.set_xlabel("Node number", fontsize=25)
        ax.legend(loc='upper right', fontsize='x-large')
        plt.show()


        # print 'point', m,'new_theta',involved_points[m], 'theta', result[m], 'time', result_count[m]
    # for n in xrange(len(result_step_by_step)):
    #     for m in result_step_by_step[n]:
    #         for k in xrange(len(result_step_by_step)):
    #             result_step_by_step[k][m]


    return result
