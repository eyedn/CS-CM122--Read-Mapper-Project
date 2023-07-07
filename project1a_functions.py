###############################################################################
#   Aydin Karatas
#   CS 122 S23
#   Project 1a
#   project1a_functions.py
###############################################################################
import project1a_classes as cl

# given a file of reads, create a list of read objects
def create_reads_dict(file, subdiv_size, min_div_size):
    reads_dict = {}
    with open(file, "r") as f:
        # define prev_read outside of the loop so we can assign pairs
        prev_read_name = ""
        for line in f:
            curr_line = line.strip()
            if '>' in curr_line:
                # names will only have read number and pair number
                read_name = curr_line[6:]
                reads_dict[read_name] = None
                continue
            else:
                # associate this sequence with the current read name
                curr_read = cl.Read(curr_line, read_name, subdiv_size, min_div_size)
                # if pair reads are detected...
                if '/' in read_name: 
                    # if pair number == 1, store the curr for future pairing
                    if read_name.split('/')[1] == '1':
                        prev_read_name = read_name
                    # pair number == 2, assign curr read and prev read as pairs
                    else:
                        curr_read.assign_pair(reads_dict[prev_read_name])
                        reads_dict[prev_read_name].assign_pair(curr_read)
                reads_dict[read_name] = curr_read
    return reads_dict

# map all reads to the genome by first mapping all subdivions of reads
def map_reads(genome, reads_dict, max_errors):
    # map the divisions of each read to the genome to produce indeces of interest
    for read in reads_dict:
        reads_dict[read].map_subdivisions(genome)
    # once the divisions have been mapped, map the reads themselves
    for read in reads_dict:
        reads_dict[read].map_read(genome.sequence, max_errors)

# find the hamming distance between two strings
def calc_hamming_dist(str_1, str_2):
    h_dist = 0
    min_len = min(len(str_1), len(str_2))
    for i in range(min_len):
        # hamming dist increases when the strings are not the same at pos i
        if str_1[i] != str_2[i]:
            h_dist += 1
    return h_dist

# once all reads are mapped, find the mutations that are likely real mutations
def concensus_mutations(reads_dict, min_occ):
    all_mutations = {}
    consensus = []
    # first, conglomerate the occurances of all mutations and errors
    for read in reads_dict:
        if reads_dict[read].mutations:
            for mut in reads_dict[read].mutations:
                try:
                    all_mutations[mut] += 1
                except:
                    all_mutations[mut] = 1
    # only consider mutations that occur >= the minimum accepted occurance lvl
    for mut in all_mutations:
        # print(f'{mut}, occ.: {all_mutations[mut]}')
        if (all_mutations[mut] >= min_occ) and (mut != None):
            consensus.append(mut)
    return consensus
    all_mutations = {}
    consensus = []
    # first, conglomerate the occurances of all mutations and errors
    for read in reads_dict:
        if reads_dict[read].mutations:
            for mut in reads_dict[read].mutations:
                try:
                    all_mutations[mut] += 1
                except:
                    all_mutations[mut] = 1
    # only consider mutations that occur >= the minimum accepted occurance lvl
    for mut in all_mutations:
        # print(f'{mut}, occ.: {all_mutations[mut]}')
        if (all_mutations[mut] >= min_occ) and (mut != None):
            consensus.append(mut)
    return consensus