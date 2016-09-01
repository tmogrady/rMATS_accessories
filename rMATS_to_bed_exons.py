#takes information from rMATS SE output files to produce a bed file of exons

import sys

file = open(sys.argv[1], 'r')

def name_stub(filename): #gets the name of the file not including the file extension
    for i in range(len(filename)-1,-1,-1):
        if filename[i] == '.':
            return filename[0:i]

name_stub = name_stub(sys.argv[1]) #gets the name of the file not including the file extension
output = open("%s_exons.bed" % name_stub, 'w') #open output file

for line in file:
    line = line.rstrip("\n") #removes carriage return
    fields = line.split("\t")
    chr = fields[3][3:]
    if fields[0] != 'ID': #skips header lines
        name = fields[1] + "|" + chr + ":" + fields[5] + "-" + fields[6]
        output.write(chr + "\t" + fields[5] + "\t" + fields[6] + "\t" + name + "|" + "up" + "\t" + fields[22] + "\t" + fields[4] + "\n")

file.close()
output.close()

#should add colour-coding for included & excluded exons