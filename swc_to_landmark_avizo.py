import os
neuron_path="SO_swc/"
record_path="FC_landmarks/"
preneuron_list=os.listdir(neuron_path)
if not os.path.isdir(record_path):
    os.mkdir(record_path)
neuron_list=[]
for i in preneuron_list:
    if i.find("SO")!=-1:
        neuron_list.append(i)
print len(neuron_list)
# print neuron_list
for neuron in neuron_list:
    if os.path.isfile(record_path+neuron[3:]+"_avizo.txt"):
        print "exist",neuron, neuron_list.index(neuron)
        continue
    print neuron, neuron_list.index(neuron)
    xyz=""
    length=""
    with open(neuron_path+neuron,"rt")as ff:
        for lines in ff:
            lines=lines[:-1]
            grouping=lines.split(" ")
            xyz=xyz+str(grouping[1])+" "+str(grouping[2])+" "+str(grouping[3])+"\n"
            length=grouping[0]
    record_lines=""
    record_lines="# AmiraMesh 3D ASCII 2.0\n"+"define Markers "+length+"\nParameters {\n"+"    ContentType \"LandmarkSet\",\n"+"    NumSets 1\n"+"}\n"+"Markers { float[3] Coordinates } @1\n"+"# Data section follows\n"+"@1\n"+xyz
    with open(record_path+neuron[3:]+"_avizo.txt","wt")as ff:
        ff.write(record_lines)




