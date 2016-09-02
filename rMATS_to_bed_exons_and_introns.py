#takes information from rMATS SE output files to produce a bed file exons with flanking introns, e.g. for motif analysis. Default is to return full introns: if the user includes a number before the file name, the script will return that number of intronic bases closest to the exon boundary (or the full intron, whichever is shorter)
#returns only one downstream and one upstream intron per exon: the longest one
#for duplicate exons (i.e. in multiple events or with multiple gene names), only one gene name and score is returned
#future improvement: color-code by score: positive, negative, or both
#USAGE: rMATS_to_bed_introns.py <number_of_bases-OPTIONAL> <rMATS_SE_file>

import sys

def name_stub(filename): #gets the name of the file not including the file extension
    for i in range(len(filename)-1,-1,-1):
        if filename[i] == '.':
            return filename[0:i]

if len(sys.argv) < 2:
    print "Not enough arguments. Usage: rMATS_to_bed_introns.py <length> <rMATS_SE_file>\n"
elif len(sys.argv) == 2:
    size = None
    file = sys.argv[1]
elif len(sys.argv) == 3:
    size = int(sys.argv[1])
    file = sys.argv[2]
else:
    print "Too many arguments. Usage: rMATS_to_bed_introns.py <length> <rMATS_SE_file>\n"

open_file = open(file, 'r')

name_stub = name_stub(file) #gets the name of the file not including the file extension

if size:
    output = open("%s_exons_and_introns_%s.bed" % (name_stub, size), 'w') #open output file
else:
    output = open("%s_exons_and_introns.bed" % name_stub, 'w') #open output file

#dictionary used for duplicate removal...
SEQ = {}

for line in open_file:
    line = line.rstrip("\n") #removes carriage return
    fields = line.split("\t")
    chr = fields[3][3:]
    if fields[0] == 'ID': #skips header lines
        continue
    name = chr + fields[4] + fields[5] + "-" + fields[6]
    if size:
        upper = int(fields[5]) - size #gets coordinates for x bp upstream...
        lower = int(fields[6]) + size #...and downstream
# to get upstream coordinate:
    if size and (upper >= int(fields[8])):  # or the whole intron, whichever is shorter.
        chrStart = upper
        if name in SEQ: #if the exon is in the dictionary already (i.e., is a duplicate)
            if int(SEQ[name][1]) > upper: # if the upstream intron in the dictionary is shorter than the current one...
                SEQ[name][1] = str(upper) #...replace it with the current, longer one, keeping other info the same
            else:
                pass
        else: # if the current exon is not in the dictionary, add it
            SEQ[name] = [chr, str(chrStart), fields[6], fields[1]+"|"+name, fields[22], fields[4]]
    else: #if the size is not specified, or is larger than the intron (i.e., full intron is being returned)
        chrStart = fields[8]
        if name in SEQ:
            if int(SEQ[name][1]) > chrStart: # if the upstream intron in the dictionary is shorter than the current one...
                SEQ[name][1] = str(chrStart) #...replace it with the current, longer one, keeping other info the same
            else:
                pass
        else: # if the current exon is not in the dictionary, add it
            SEQ[name] = [chr, str(chrStart), fields[6], fields[1]+"|"+name, fields[22], fields[4]]
# to get downstream coordinate:
    if size and (lower <= int(fields[9])):
        chrEnd = lower
        if int(SEQ[name][2]) < chrEnd: #don't need "if name in SEQ" because name will always be in SEQ from above
            SEQ[name][2] = str(chrEnd)
        else:
            pass
    else:
        chrEnd = fields[9]
        if int(SEQ[name][2]) < chrEnd:
            SEQ[name][2] = str(chrEnd)
        else:
            pass

for key in SEQ:
    output.write("\t".join(SEQ[key]) + "\n")

open_file.close()
output.close()