# to get population number to determine proportions: count lines? do we want enrichment per intron or per exon? For introns, there should be 2 sequences per exon. For UTRs there is not a consistent number of UTRs per gene, becuase of multiple transcripts. For now maybe keep total as user input
#but, if we're not using number of sequences as the total population number, need to adjust how we count. Right now if an exon has the same motif upstream and downstream, it's counted twice. Maybe find a way to concatenate sequences corresponding to the same exon (or gene, or whatever) without getting motifs from the connection point?

#USAGE: <script.py> <desired_motif_length> <sample_file> <control_file>

import sys,itertools,operator
import matplotlib.pyplot as plt

l = int(sys.argv[1]) #length of motifs to investigate

def generate_motifs(): #to generate all possible motifs of the given length
    MOTIFS = {}
    for item in itertools.product('ATCG', repeat=l):
        MOTIFS[''.join(item)] = 0
    return MOTIFS

def get_motifs(seq):
    return list(set([seq[i:i+l] for i in range(0,len(seq)-l+1)])) # to find all motifs in provided sequences. To improve: use a more efficient way of doing this than generating the whole list and then removing duplicates

def count_motif_seqs(file, DICTIONARY): # to get counts of how many sequences have each motif
    line_count = 0
    for line in file:
        line = line.rstrip('\n')
        if line[0] == ">":
            line_count += 1
        else:
            seq_motifs = get_motifs(line) #finds motifs in the sequence line
            for motif in seq_motifs:
                if motif.count('N') == 0: #change this so that motifs must have A C T G only
                    DICTIONARY[motif] += 1 #adds this sequence to the dictionary count for the detected sequences
            #print seq_motifs. DO I need to return this?
    return line_count

def z_score(DICT1, DICT2, n1, n2):
    Z_DICT = {}
    n1 = float(n1)
    n2 = float(n2)
    for key in DICT1:
        if DICT1[key] != 0 or DICT2[key] != 0:
            p1 = DICT1[key] / (n1)
            p2 = DICT2[key] / (n2)
            p = (DICT1[key]+DICT2[key]) / (n1+n2)
            z = (p1-p2) / (p*(1-p)*(1/n1 + 1/n2))**0.5
        else:
            z = 0.0
        Z_DICT[key] = z
    return Z_DICT

#call the function to generate the list of possible motifs
SAMPLE_MOTIFS = generate_motifs()
CONTROL_MOTIFS = generate_motifs()
#print len(SAMPLE_MOTIFS) #to make sure it's working. E.g. there are 1024 possible pentamers

#call the function to find motifs present in each sequence, from both the sample file and the control file, and count them
print "getting sample motifs"
sample_file = open(sys.argv[2], 'r')
sample_seq_number = count_motif_seqs(sample_file, SAMPLE_MOTIFS)
#print sample_seq_number
sample_file.close()

print "getting control motifs"
control_file = open(sys.argv[3], 'r')
control_seq_number = count_motif_seqs(control_file, CONTROL_MOTIFS)
#print control_seq_number
control_file.close()

#print "SAMPLE"
#for key in SAMPLE_MOTIFS:
#    print key + " " + str(SAMPLE_MOTIFS[key])
#
#print "CONTROL"
#for key in CONTROL_MOTIFS:
#    print key + " " + str(CONTROL_MOTIFS[key])

#calculate the z-score for enrichment of each motif
print "calculating z-scores"
Z_DICT = z_score(SAMPLE_MOTIFS, CONTROL_MOTIFS, sample_seq_number, control_seq_number)
#print Z_DICT

#sort the z-scores to list motifs in descending order
sorted_Z_DICT = sorted(Z_DICT.items(), key=operator.itemgetter(1), reverse=True)
#print sorted_Z_DICT

output = open("%dmer_z-scores.txt" % l, 'w') #add a header line with names of script and input files
output.write("%dmer z-scores from motif_enrichment_finder_per_sequence.py\ninput files: %s, %s\n" % (l,sys.argv[2],sys.argv[3]))
output.write('\n'.join('%s %s' % x for x in sorted_Z_DICT))
output.close

#plot the z-scores as a histogram. Maybe should make this a function and call it?
z_to_plot = Z_DICT.values()
#print z_to_plot
plt.hist(z_to_plot, bins=25, color='#483D8B')
plt.title("%dmer Z-score distribution" % l)
plt.xlabel("%dmer Z-score" % l)
plt.ylabel("Number of sequences with %dmer" % l)
plt.savefig("%dmer_histogram.png" % l)
#
#plt.hist(z_to_plot, bins=25, normed=True, color='#008B8B')
#plt.title("%dmer Z-score distribution" % l)
#plt.xlabel("%dmer Z-score" % l)
#plt.ylabel("Probability of a sequence with %dmer" % l)
#plt.savefig("%dmer_normed_histogram.png" % l)
