import os
import re
# url=[]
file_list=[]
path="xbrain_new_20000/warpResult/"
path_record="xbrain_new_20000_FC/"
for a,b,c in os.walk(path):
    for i in c:
        if i.find("FC")!=-1:
           file_list.append(i)
try:
    os.mkdir(path_record)
except:
    pass
for file in file_list:
    landmarks=[]
    with open (path+file,"rt") as ff:
        check=0
        for line in ff:
            if check==0:
                if line[0]=="@" and line[1]=="1":
                    check=1
            elif check==1 and line[0]=="@" and line[1]=="2":
                break
            else:
                line=line[:-1]
                grouping=line.split(" ")
                if len(grouping)!=3:
                    print grouping
                    break;
                landmarks.append([float(grouping[0]),float(grouping[1]),float(grouping[2])])
        print file
        # print landmarks
    with open(path_record+file,"wt")as fp:
        for [x,y,z] in landmarks:
            fp.writelines(str(x)+" "+str(y)+" "+str(z)+"\n")
