#takes information from rMATS SE output files to produce a bed file of introns flanking the exons. User can indicate whether they want the full intron or only the specified number of basepairs nearest the exon.

import sys

size = int(sys.argv[1]) # fix this so that default is full intron
file = open(sys.argv[2], 'r')

def name_stub(filename): #gets the name of the file not including the file extension
    for i in range(len(filename)-1,-1,-1):
        if filename[i] == '.':
            return filename[0:i]

name_stub = name_stub(sys.argv[2]) #gets the name of the file not including the file extension
output = open("%s_introns_%s.bed" % (name_stub, size), 'w') #open output file

for line in file:
    line = line.rstrip("\n") #removes carriage return
    fields = line.split("\t")
    chr = fields[3][3:]
    if fields[0] == 'ID': #skips header lines
        pass
    name = fields[1] + "|" + chr + ":" + fields[5] + "-" + fields[6] #later, add in "ex" and "in"
    if fields[4] == "+":
        upper = int(fields[5]) - size #gets coordinates for 200 bp upstream...
        if upper >= int(fields[8]):  # or the whole intron, whichever is shorter. Is this int() necessary? Does fields[8] start as a string?
            output.write(chr + "\t" + str(upper) + "\t" + fields[5] + "\t" + name + "|" + "up" + "\t" + fields[22] + "\t" + fields[4] + "\n")
        else:
            output.write(chr + "\t" + fields[8] + "\t" + fields[5] + "\t" + name + "|" + "up" + "\t" + fields[22] + "\t" + fields[4] + "\n")
        lower = int(fields[6]) + size
        if lower <= int(fields[9]):
            output.write(chr + "\t" + fields[6] + "\t" + str(lower) + "\t" + name + "|" + "down" + "\t" + fields[22] + "\t" + fields[4] + "\n")
        else:
            output.write(chr + "\t" + fields[6] + "\t" + fields[9] + "\t" + name + "|" "down" + "\t" + fields[22] + "\t" + fields[4] + "\n")
    elif fields[4] == "-":
        upper = int(fields[6]) + size
        name = fields[1] + "|" + chr + ":" + fields[5] + "-" + fields[6]
        if upper <= int(fields[9]): #not sure if int is necessary
            output.write(chr + "\t" + fields[6] + "\t" + str(upper) + "\t" + name + "|" + "up" + "\t" + fields[22] + "\t" + fields[4] + "\n")
        else:
            output.write(chr + "\t" + fields[6] + "\t" + fields[9] + "\t" + name + "|" + "up" + "\t" + fields[22] + "\t" + fields[4] + "\n")
        lower = int(fields[5]) - size
        if lower >= int(fields[7]):
            output.write(chr + "\t" + str(lower) + "\t" + fields[5] + "\t" + name + "|" + "down" + "\t" + fields[22] + "\t" + fields[4] + "\n")
        else:
            output.write(chr + "\t" + fields[8] + "\t" + fields[5] + "\t" + name + "|" + "down" + "\t" + fields[22] + "\t" + fields[4] + "\n")

#need to:
#figure out 0- and 1-based issues for starts and ends



file.close()
output.close()