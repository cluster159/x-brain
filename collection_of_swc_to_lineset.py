import os
path_similar_neuron="result_xbrain2FC_neuron_similariy0p4_distance6_sig2um/"
neuron_path="SO_swc/"
record_path="result_rebuild_0p2_xbrain2FC_neuron_similariy0p4_distance6_sig2um_lineset/"
brain_list=os.listdir(path_similar_neuron)
preneuron_list=os.listdir(neuron_path)

if not os.path.isdir(record_path):
    os.mkdir(record_path)
for brain in brain_list:
    neuron_list=[]
    with open(path_similar_neuron+brain,"rt")as ff:
        for line in ff:
            line=line[:-1]
            grouping=line.split(" ")
            print grouping
            similarity=float(grouping[2])
            neuron_name=grouping[0]
            if similarity>=0.2:
                neuron_list.append(neuron_name)
            else:
                break
    print len(neuron_list)
    print "search over", brain_list.index(brain)
    length = 0
    line_number = 0
    nodes = []
    lines = []
    line = []
    for neuron in neuron_list:
        print neuron, neuron_list.index(neuron)
        with open(neuron_path +"SO_"+ neuron, "rt")as ff:
            for data in ff:
                data = data[:-1]
                grouping = data.split(" ")
                nodes.append([grouping[1], grouping[2], grouping[3]])
                if int(grouping[4]) != -1 and (int(grouping[0]) - int(grouping[4])) != 1:
                    line.append(-1)
                    lines.append(line)
                    line = []
                    if grouping[1] != nodes[int(grouping[4]) - 1][0] and grouping[2] != nodes[int(grouping[4]) - 1][
                        1] and grouping[3] != nodes[int(grouping[4]) - 1][2]:
                        line.append(int(grouping[4]) + length - 1)
                line.append(int(grouping[0]) +length - 1)
        line.append(-1)
        lines.append(line)
        length = int(grouping[0])+length
        line_number = len(lines)+line_number
        line = []
    record_lines = ""
    record_lines = "# AmiraMesh 3D ASCII 2.0\n" + "define Lines " + str(line_number*10) + "\nnVertices " + str(
            length) + "\nParameters {\n" + "    ContentType \"HxLineSet\"\n}\n" + "Lines { int LineIdx } @2\n" + "Vertices { float[3] Coordinates } @1\n" + "# Data section follows\n" + "@1\n"
    with open(record_path + brain + "_FC_avizo.txt", "wt")as ff:
        ff.write(record_lines)
        for node in nodes:
            ff.writelines(str(node[0]) + " " + str(node[1]) + " " + str(node[2]) + "\n")
        ff.writelines("@2\n")
        for line in lines:
            for node in line:
                ff.writelines(str(node) + " ")
            ff.writelines("\n")
    del (line)
    del (lines)
    del (record_lines)





