# Given a bed file, finds overlapping features and randomly selects one to return. Bed file must be sorted by chromosome then start and end position
# to improve: make script sort the file first in case user hasn't
# USAGE: overlapping_bed_remover.py <bed_file.bed>

import sys, random

def name_stub(filename): #gets the name of the file not including the file extension
    for i in range(len(filename)-1,-1,-1):
        if filename[i] == '.':
            return filename[0:i]

file = sys.argv[1]

open_file = open(file, 'r')
name_stub = name_stub(file)

output = open("%s_no_overlapping.bed" % name_stub, 'w')

cluster = []

for line in open_file:
    line = line.rstrip("\n")
    fields = line.split("\t")
    if len(cluster) == 0: #adds values for first feature to the cluster
        cluster.append([fields[0], fields[5], fields[1], fields[2]])
    else: #checks to see if current feature overlaps other features in the cluster
        chrStart = min([cluster[x][2] for x in range(0, len(cluster))])
        chrEnd = max([cluster[x][3] for x in range(0, len(cluster))])
        if fields[0]==cluster[0][0] and fields[5]== cluster[0][1] and ((chrStart <= fields[1] <= chrEnd) or (chrStart <= fields[2] <= chrEnd)):
            cluster.append([fields[0], fields[5], fields[1], fields[2]])
        else: #randomly choose one of the features in the cluster to print out
            random_index = random.randint(0, len(cluster)-1)
            bed_line = cluster[random_index][0] + "\t" + cluster[random_index][2] + "\t" + cluster[random_index][3] + "\t" + cluster[random_index][0]+":"+cluster[random_index][2]+"-"+cluster[random_index][3] + "\t1\t" + cluster[random_index][1] + "\n"
            output.write(bed_line)
            cluster = []

random_index = random.randint(0, len(cluster)-1) #print out last feature
bed_line = cluster[random_index][0] + "\t" + cluster[random_index][2] + "\t" + cluster[random_index][3] + "\t" + cluster[random_index][0]+":"+cluster[random_index][2]+"-"+cluster[random_index][3] + "\t1\t" + cluster[random_index][1] + "\n"

open_file.close()
output.close()