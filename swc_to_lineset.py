import os
neuron_path="SO_swc/"
record_path="FC_linesets/"
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
    line_number=0
    # if os.path.isfile(record_path+neuron[3:]+"_avizo.txt"):
    #     print "exist",neuron, neuron_list.index(neuron)
    #     continue
    print neuron, neuron_list.index(neuron)
    nodes=[]
    length=0
    lines=[]
    line=[]
    with open(neuron_path+neuron,"rt")as ff:
        for data in ff:
            data=data[:-1]
            grouping=data.split(" ")
            nodes.append([grouping[1],grouping[2],grouping[3]])
            if int(grouping[4])!=-1 and (int(grouping[0])-int(grouping[4]))!=1:
                line.append(-1)
                lines.append(line)
                # print line
                line=[]
                if grouping[1]!=nodes[int(grouping[4])-1][0] and grouping[2]!=nodes[int(grouping[4])-1][1] and grouping[3]!=nodes[int(grouping[4])-1][2]:
                    line.append(int(grouping[4])-1)
            line.append(int(grouping[0])-1)
    line.append(-1)
    lines.append(line)
    length=grouping[0]
    line_number=len(lines)+int(length)+500
    record_lines=""
    record_lines="# AmiraMesh 3D ASCII 2.0\n"+"define Lines "+str(line_number)+"\nnVertices "+str(length)+"\nParameters {\n"+"    ContentType \"HxLineSet\"\n}\n"+"Lines { int LineIdx } @2\n"+"Vertices { float[3] Coordinates } @1\n"+"# Data section follows\n"+"@1\n"
    with open(record_path+neuron[3:]+"_avizo.txt","wt")as ff:
        ff.write(record_lines)
        for node in nodes:
            ff.writelines(str(node[0])+" "+str(node[1])+" "+str(node[2])+"\n")
        ff.writelines("@2\n")
        for line in lines:
            for node in line:
                ff.writelines(str(node)+" ")
            ff.writelines("\n")





