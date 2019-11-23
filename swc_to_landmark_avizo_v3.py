import os
neuron_path="20191115_new_bin1_xbrain/19000/"
record_path="20191115_new_bin1_xbrain_19000_xyz/"
preneuron_list=os.listdir(neuron_path)
print preneuron_list
if not os.path.isdir(record_path):
    os.mkdir(record_path)
neuron_list=[]
for i in preneuron_list:
    if i.find("Results")!=-1:
        neuron_list.append(i[:-8]+".swc")
print len(neuron_list)
# print neuron_list
for neuron in neuron_list:
    if os.path.isfile(record_path+neuron[:-4]+"_xyz.txt"):
        print "exist",neuron, neuron_list.index(neuron)
        continue
    with open(record_path+neuron[:-4]+"_xyz.txt","wt") as ff:
        ff.write("")
    print neuron, neuron_list.index(neuron)
    x=[]
    y=[]
    z=[]
    length=""
    with open(neuron_path+neuron[:-4]+"-Results/"+neuron,"rt")as ff:
        for lines in ff:
            lines=lines[:-1]
            grouping=lines.split(" ")
            x.append(grouping[2])
            y.append(grouping[3])
            z.append(grouping[4])
            length=grouping[0]
    length=str(int(grouping[0]*10))
    record_lines = ""
    record_lines = "# AmiraMesh 3D ASCII 2.0\n" + "define Markers " + length + "\nParameters {\n" + "    ContentType \"LandmarkSet\",\n" + "    NumSets 1\n" + "}\n" + "Markers { float[3] Coordinates } @1\n" + "# Data section follows\n" + "@1\n"
    with open(record_path+neuron[:-4]+"_xyz.txt","wt")as ff:
        ff.write(record_lines)
        for i in xrange(len(x)):
            ff.writelines(x[i]+" "+y[i]+" "+z[i]+"\n")
            




