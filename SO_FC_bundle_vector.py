#! /usr/bin/env python
# coding: utf-8
import numpy as np
import math as math
import time
from collections import OrderedDict

def pos_to_int(x): #convert to int(x*10)
    return int(round(10*x))
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

neuro_space=[]
neuro_space_list=[]


def name_to_index(name): #convert name to int
    return nn.namelist.index(name)


def detect_neuron(inputx, inputy, inputz,input_vec,detect_list,brain_dictionary,cube_size): #version for input is a single point; datapoints is an OrderedDict
    involved_neurons = OrderedDict() #record neuron with its possibility ((key:n名稱):(value: p機率),(n):(p),...)
    for i in detect_list:
        x=i[0]+inputx
        y=i[1]+inputy
        z=i[2]+inputz
        p=i[3]
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
        if (x_index,y_index,z_index) in brain_dictionary:
            if (x,y,z) in brain_dictionary[(x_index,y_index,z_index)]:
                for k in brain_dictionary[(x_index,y_index,z_index)][(x,y,z)]:
                    # print k
                    dot = float(k[3] * input_vec[0] + k[4] * input_vec[1] + k[5] * input_vec[2])
                    # dot = float(dot / 100)
                    # print "HERE"
                    # print input_vec
                    # print k
                    # print dot
                    dot = np.clip(np.abs(dot), 0, 1)
                    if ((k[0],k[1])) not in involved_neurons:
                        involved_neurons[(k[0],k[1])] = dot
                    else:
                        if involved_neurons[(k[0],k[1])] < dot:
                            involved_neurons[(k[0],k[1])] = dot
    return involved_neurons


def inference(inputpoints,detect_list,brain_dictionary,cube_size): #assume inputpoints i s a list of lists/tuples; datapoints is an OrderedDict，這個沒包含vector
    # print "infer"
    global result
    result={}
    result_count={}
    neuron_number = 0
    count=0
    for i in xrange(len(inputpoints)):
        input_vec=(inputpoints[i][6],inputpoints[i][7],inputpoints[i][8])
        involved_points = detect_neuron(inputx=inputpoints[i][0],inputy=inputpoints[i][1],inputz=inputpoints[i][2],input_vec=input_vec,detect_list=detect_list,brain_dictionary=brain_dictionary,cube_size=cube_size)
        neuron_new_number=len(involved_points)
        if neuron_new_number == 0 and neuron_number ==0:
            continue
        else:
            if neuron_new_number>neuron_number:
                neuron_number = neuron_new_number
            count = count + 1
            if count == 1:
                result = involved_points
                for m in result:
                    result_count[m]=1
            else:
                for m in involved_points:
                    if m in result:
                        # print 'm=',m
                        result[m]=result[m]+involved_points[m]
                        result_count[m] = result_count[m] + 1

                    else:
                        result[m] = involved_points[m]
                        result_count[m] =1

    #             四種狀況，在result(舊的)和involved(新的)都出現，只出現在result，只出現在involved，兩個都沒出現過
    #             都出現的不管是哪個人跑迴圈都可以偵測到，只出現在result的，不需做處理，因為直接加上0，而只出現在inovlved，則需要加上新的值(過去的都自動設為0，所以都不用處理)，如果是兩個都沒出現，當然不用處理
    #              →因此只需要跑involved的迴圈就可!!!!!!!!!!
    # print 'number',neuron_new_number,neuron_number
    # print 'all',result
    # del neuro_space_list
    # del neuro_space
    neuron_length=len(inputpoints)
    for i in result:
        result[i]=result[i]/result_count[i]
    return result

