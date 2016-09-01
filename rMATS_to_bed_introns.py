#takes information from rMATS SE output files to produce a bed file of introns flanking the exons. Default is to return the full intron: if the user includes a number before the file name, the script will return that number of intronic bases closest to the exon boundary (or the full intron, whichever is shorter)
#returns only one downstream and one upstream intron per exon: the longest one
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
    output = open("%s_introns_%s.bed" % (name_stub, size), 'w') #open output file
else:
    output = open("%s_introns.bed" % name_stub, 'w') #open output file

#dictionaries used for duplicate removal...
UPSTREAM = {}
DOWNSTREAM = {}
#but! Can nest dictionaries. So try that later to make things more concise

for line in open_file:
    line = line.rstrip("\n") #removes carriage return
    fields = line.split("\t")
    chr = fields[3][3:]
    if fields[0] == 'ID': #skips header lines
        pass
    name = chr + fields[4] + fields[5] + "-" + fields[6] #need to change this to not include the name  #later, add in "ex" and "in"
    if fields[4] == "+":
        if size:
            upper = int(fields[5]) - size #gets coordinates for x bp upstream...
            lower = int(fields[6]) + size #...and downstream
        if size and (upper >= int(fields[8])):  # or the whole intron, whichever is shorter. Is this int() necessary? Does fields[8] start as a string?
            intron = [chr, upper, fields[5], fields[1]+"|"+name+"_up", fields[22], fields[4]]
            if name in UPSTREAM:#checks to see if this exon has been encountered before
                if int(UPSTREAM[name][2]) - int(UPSTREAM[name][1]) < int(intron[2]) - int(intron[1]): #if so, checks to see if the current intron is longer than the previously see intron
                    UPSTREAM[name] = intron #if so, replaces the previous intron with the current intron
                else:
                    pass #ignores current intron if it's shorter than one previously seen for that exon
            else:
                UPSTREAM[name] = intron #if this exon hasn't been seen before, adds it to the dictionary with its intron
#output.write(chr + "\t" + str(upper) + "\t" + fields[5] + "\t" + name + "|" + "up" + "\t" + fields[22] + "\t" + fields[4] + "\n")
        else:
            intron = [chr, fields[8], fields[5], fields[1]+"|"+name+"_up", fields[22], fields[4]]
            if name in UPSTREAM:#checks to see if this exon has been encountered before
                if int(UPSTREAM[name][2]) - int(UPSTREAM[name][1]) < int(intron[2]) - int(intron[1]): #if so, checks to see if the current intron is longer than the previously see intron
                    UPSTREAM[name] = intron #if so, replaces the previous intron with the current intron
                else:
                    pass #ignores current intron if it's shorter than one previously seen for that exon
            else:
                UPSTREAM[name] = intron #if this exon hasn't been seen before, adds it to the dictionary with its intron
#output.write(chr + "\t" + fields[8] + "\t" + fields[5] + "\t" + name + "|" + "up" + "\t" + fields[22] + "\t" + fields[4] + "\n")
        if size and (lower <= int(fields[9])):
            intron = [chr, fields[6], lower, fields[1]+"|"+name+"_down", fields[22], fields[4]]
            if name in DOWNSTREAM:
                if int(DOWNSTREAM[name][2]) - int(DOWNSTREAM[name][1]) < int(intron[2]) - int(intron[1]):
                    DOWNSTREAM[name] = intron
                else:
                    pass
            else:
                DOWNSTREAM[name] = intron
#output.write(chr + "\t" + fields[6] + "\t" + str(lower) + "\t" + name + "|" + "down" + "\t" + fields[22] + "\t" + fields[4] + "\n")
        else:
            intron = [chr, fields[6], fields[9], fields[1]+"|"+name+"_down", fields[22], fields[4]]
            if name in DOWNSTREAM:
                if int(DOWNSTREAM[name][2]) - int(DOWNSTREAM[name][1]) < int(intron[2]) - int(intron[1]):
                    DOWNSTREAM[name] = intron
                else:
                    pass
            else:
                DOWNSTREAM[name] = intron
            #output.write(chr + "\t" + fields[6] + "\t" + fields[9] + "\t" + name + "|" "down" + "\t" + fields[22] + "\t" + fields[4] + "\n")
    elif fields[4] == "-":
        if size:
            upper = int(fields[6]) + size
            lower = int(fields[5]) - size
        #name = fields[1] + "|" + chr + ":" + fields[5] + "-" + fields[6]
        if size and (upper <= int(fields[9])): #not sure if int is necessary
            intron = [chr, fields[6], upper, fields[1]+"|"+name+"_up", fields[22], fields[4]]
            if name in UPSTREAM:
                if int(UPSTREAM[name][2]) - int(UPSTREAM[name][1]) < int(intron[2]) - int(intron[1]): #if so, checks to see if the current intron is longer than the previously see intron
                    UPSTREAM[name] = intron #if so, replaces the previous intron with the current intron
                else:
                    pass #ignores current intron if it's shorter than one previously seen for that exon
            else:
                UPSTREAM[name] = intron #if this exon hasn't been seen before, adds it to the dictionary with its intron
        #output.write(chr + "\t" + fields[6] + "\t" + str(upper) + "\t" + name + "|" + "up" + "\t" + fields[22] + "\t" + fields[4] + "\n")
        else:
            intron = [chr, fields[6], fields[9], fields[1]+"|"+name+"_up", fields[22], fields[4]]
            if name in UPSTREAM:
                if int(UPSTREAM[name][2]) - int(UPSTREAM[name][1]) < int(intron[2]) - int(intron[1]): #if so, checks to see if the current intron is longer than the previously see intron
                    UPSTREAM[name] = intron #if so, replaces the previous intron with the current intron
                else:
                    pass #ignores current intron if it's shorter than one previously seen for that exon
            else:
                UPSTREAM[name] = intron #if this exon hasn't been seen before, adds it to the dictionary with its intron
                    #output.write(chr + "\t" + fields[6] + "\t" + fields[9] + "\t" + name + "|" + "up" + "\t" + fields[22] + "\t" + fields[4] + "\n")
        if size and (lower >= int(fields[7])):
            intron = [chr, lower, fields[5], fields[1]+"|"+name+"_down", fields[22], fields[4]]
            if name in DOWNSTREAM:
                if int(DOWNSTREAM[name][2]) - int(DOWNSTREAM[name][1]) < int(intron[2]) - int(intron[1]):
                    DOWNSTREAM[name] = intron
                else:
                    pass
            else:
                DOWNSTREAM[name] = intron
        #output.write(chr + "\t" + str(lower) + "\t" + fields[5] + "\t" + name + "|" + "down" + "\t" + fields[22] + "\t" + fields[4] + "\n")
        else:
            intron = [chr, fields[8], fields[5], fields[1]+"|"+name+"_down", fields[22], fields[4]]
            if name in DOWNSTREAM:
                if int(DOWNSTREAM[name][2]) - int(DOWNSTREAM[name][1]) < int(intron[2]) - int(intron[1]):
                    DOWNSTREAM[name] = intron
                else:
                    pass
            else:
                DOWNSTREAM[name] = intron
#output.write(chr + "\t" + fields[8] + "\t" + fields[5] + "\t" + name + "|" + "down" + "\t" + fields[22] + "\t" + fields[4] + "\n")

#need to:
#figure out how to deal with duplicates

#output.write('\n'.join('%s %s' % x for x in sorted_Z_DICT))

for key in UPSTREAM:
    output.write("\t".join(UPSTREAM[key]) + "\n")
for key in DOWNSTREAM:
    output.write("\t".join(DOWNSTREAM[key]) + "\n")

open_file.close()
output.close()