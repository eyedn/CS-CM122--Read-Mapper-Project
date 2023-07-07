###############################################################################
#   Aydin Karatas
#   CS 122 S23
#   Project 1a
#   project1a_classes.py
###############################################################################
import project1a_functions as fx

# a genome stores the whole sequence and a burrows wheeler transform of itself
class Genome:
    def __init__(self, file):
        self.sequence = ""
        self.bwt = None
        # generate sequence
        self.read_file(file)
        # create bwt transform of sequence; creation of bwt uses a suffix array
        self.bwt_sequence()
    
    # read file file that contains the geneom
    def read_file(self, file):
        with open(file, "r") as f:
            for line in f:
                curr_line = line.strip()
                if '>' in curr_line:
                    continue
                self.sequence += curr_line.upper()

    # after the sequence has been established, create the bwt of the sequence
    def bwt_sequence(self):
        # add sentinel character $ to end of sequence
        self.sequence += '$'
        self.bwt = BWT(self.sequence)

        # given a BWT and a pattern, return the indeces where the pattern is found
    
    # map a pattern to the BWT object and return the top and bottom indeces
    def bw_matching(self, pattern):
        # start with top and bottom indeces at the start and end of the sequence
        top = 0
        bottom = len(self.bwt.last) - 1
        while top <= bottom:
            if pattern:
                # isolate the last char in the pattern
                symbol = pattern[-1]
                pattern = pattern[:-1]
                # if the symbol is found...
                if symbol in self.bwt.last[top:bottom+1]:
                    # use the last col to find the first and last index of symbol
                    top_idx_in_last = self.bwt.last[top:bottom+1].index(symbol) + top
                    bottom_idx_in_last = self.bwt.last[top:bottom+1].rindex(symbol) + top
                    # these indeces will be used to redefine top and bottom
                    top = self.bwt.last_idx[top_idx_in_last] 
                    bottom = self.bwt.last_idx[bottom_idx_in_last]
                    # return an invalid index if the pattern is not found
                else:
                    return None, None
            # once we have reached the end of the pattern...
            else:
                # return the top and bottom indeces; this will be used on
                # the suffix array to find the indeces where this pattern occurs
                return top, bottom

class Suffix_Array:
    def __init__(self, text):
        self.indeces = []
        self.construct(text)

    def construct(self, text):
        # Use a generator expression to generate the suffixes
        suffixes = (text[i:] for i in range(len(text)))
        # Use quicksort to sort the suffixes
        sorted_suffixes = sorted(suffixes)
        # Use slicing to get the indices of the sorted suffixes
        self.indeces = [len(text) - len(suffix) for suffix in sorted_suffixes]

# burrow's wheeler transform
class BWT:
    def __init__(self, text):
        # store the first and last cols of the BWT(T) matrix
        self.first = ''
        self.last = ''
        # store the suffix array to index where each row of the matrix starts
        self.suffix_array = []
        # generate first and last cols; generate suffix array
        self.bwt_from_suffix(text)
        # store indeces of the first col and how they relate to the last row
        self.first_idx = []
        self.last_idx = []
        # generate indeces for matching
        self.create_indeces()

    # generate a suffix array for the bwt; bwt is modelled off a suffix array
    def bwt_from_suffix(self, text):
        self.suffix_array = Suffix_Array(text)
        # last col is the letter right before the letters index in the array
        self.last = "".join(text[idx - 1] for idx in self.suffix_array.indeces)
        # first col is the last col lexicographically
        self.first = "".join(sorted(self.last))

    # once the first and last cols are established, generate i and lf-i indeces
    def create_indeces(self):
        # indeces for first are 1 through max(first length)
        self.first_idx = list(range(0, len(self.first)))
        # indeces for last require matching with first
        first_idx_count = {}
        for i,char in enumerate(self.first):
            try:
                first_idx_count[char].append(i)
            except:
                first_idx_count[char] = [i]
        for char in self.last:
            self.last_idx.append(first_idx_count[char][0])
            first_idx_count[char] = first_idx_count[char][1:]

# a read while keep track of its subdivisions and where those subdivisions map
class Read:
    def __init__(self, sequence, read_label, divisions_length, min_div_size):
        self.sequence = sequence.upper()
        self.label = read_label
        self.divisions = []
        self.division_idx = []
        self.division_length = int(divisions_length)
        self.create_divisions(min_div_size)
        # if a read is paried, a read will know who its pair is
        self.pair = None
        # once mapped, a read will keep track of its genome index and mutations
        self.mapped_idx = None
        self.mutations = []
    
    # assign anothe read object as the partner of this read
    def assign_pair(self, partner):
        self.pair = partner

    # divide the read into 
    def create_divisions(self, min_div_size):
        # divide sequence into subdivisions of length division_length
        for idx in range(0, len(self.sequence), self.division_length):
            self.divisions.append(self.sequence[idx:idx+self.division_length])
        # check if there are at least 2 subdivs the last division is too small
        if len(self.divisions) >= 2 and len(self.divisions[-1]) < min_div_size:
            # add the last subdiv (thats too small) to the 2nd-to-last subdiv
            self.divisions[-2] += self.divisions[-1]
            self.divisions.pop()

    # map all subdivisions to the genome and store their indeces
    def map_subdivisions(self, genome):
        for subdiv in self.divisions:
            # genearte the top and bottom indeces for the bwt suffix array
            top, bottom = genome.bw_matching(subdiv)
            # if the indeces were found...
            if top:
                # take the indeces from top to bottom of the suffix array
                indeces = genome.bwt.suffix_array.indeces[top:bottom+1]
                self.division_idx.append(indeces)
            else:
                self.division_idx.append([])

    # once the subdivision are mapped, use their mapped indeces to map the read
    def map_read(self, genome, max_errors):
        # for each subdivision...
        for i, subdiv_idx in enumerate(self.division_idx):
            # for each position index for this subdivision...
            for idx in subdiv_idx:
                # define the start position in genome and length of read
                # the starting position needs to be adjusted based 
                # on which subdivion this index comes form
                s = idx - i*self.division_length
                l = len(self.sequence)
                # if staring idx is 0, this subvision was mapped incorrectly
                if s < 0:
                    continue
                # if the hamming distance is <= the max allowed errors,
                # we define this index as this read's index
                h_dist = fx.calc_hamming_dist(self.sequence, genome[s:s+l])
                if h_dist <= int(max_errors):
                    self.mapped_idx = s
                    self.find_snps(genome[s:s+l], s, 0)
                    return
                # if the hamming dist is greater than the allowed errors,
                # check for an insertions and try and realign the genome
                # retry the hamming distance
                div_align = self.find_related_indeces(idx, i)
                reference, found_indel, indel_shift = self.find_indels(genome, s, i, div_align)
                h_dist = fx.calc_hamming_dist(self.sequence, reference)
                if h_dist <= int(max_errors):
                    self.mapped_idx = s
                    self.find_snps(reference, s, indel_shift)
                    if found_indel != None:
                        self.mutations.append(found_indel)
                    return

    # find associated indeces of a given index, returning how much they deviate
    def find_related_indeces(self, index, subdiv_num):
        indeces = []
        for s, subdiv in enumerate(self.division_idx):
            # return this index if we are at the subdiv of the index
            if s == subdiv_num:
                indeces.append(0)
                continue
            # return an invalid index if the subdiv is empty
            if not subdiv:
                indeces.append(None)
                continue
            # determine the ideal related index for this subdiv
            steps = abs(s - subdiv_num)
            if s < subdiv_num:
                search = index - steps*self.division_length
            else:
                search = index + steps*self.division_length
            # find the closet value so the ideal index
            closest_num = min(subdiv, key=lambda x: abs(x - search))
            indeces.append(closest_num - search)
        return indeces

    # find locations that show signs of an indel; return adjusted genome window
    def find_indels(self, genome, mapped_idx, div_idx, div_alignment):
        # if an possible indel is found, real values will be given to these
        shift = None
        is_insertion = None
        indel_div = None
        # determine whether an indel exists and in which subdivision
        for i, div_val in enumerate(div_alignment):
            # skip the current division; skip divsions with no mapping
            if i == div_idx:
                continue
            # record the first none division as the indel div
            if (div_val == None) and (not indel_div):
                indel_div = i
                continue
            # check the case where there is a shift in the indeces
            if div_val and div_val != 0:
                if (i > div_idx and div_val > 0) or (i < div_idx and div_val < 0):
                    is_insertion = False
                else:
                    is_insertion = True
                shift = abs(div_val)
            # define the divsion with the indel if not already defined
            if shift and not indel_div:
                indel_div = i
        # if no shift is detected return original genome window and no indels
        if not shift:
            return genome[mapped_idx:mapped_idx+len(self.sequence)], None, 0
        # check the division of interest for the indel
        indel_s = mapped_idx+indel_div*self.division_length
        indel_l = len(self.divisions[indel_div])
        reference = genome[indel_s:indel_s+indel_l]
        for i in range(len(self.divisions[indel_div])):
            # find mismatch in indel div to find the indel
            if self.divisions[indel_div][i] != reference[i]:
                # if insertion, add the query insertion to the reference and
                # delete the last base in the reference
                if is_insertion:
                    adjusted_genome = genome[mapped_idx:indel_s+i] + self.divisions[indel_div][i:i+shift] + genome[indel_s+i:mapped_idx+len(self.sequence)-1]
                    found_indel = f'>I{mapped_idx+indel_div*self.division_length+i-1} {self.divisions[indel_div][i:i+shift]}'
                    adjust_for_shift = 0 - shift
                # if deletion, delete the base not found in the query and
                # add the next base in the geneom not included in the reference
                else:
                    adjusted_genome = genome[mapped_idx:indel_s+i] + genome[indel_s+i+shift:mapped_idx+len(self.sequence)+1]
                    found_indel = f'>D{mapped_idx+indel_div*self.division_length+i-1} {reference[i]}'
                    adjust_for_shift = shift
                return adjusted_genome, found_indel, adjust_for_shift
                break
        return genome[mapped_idx:mapped_idx+len(self.sequence)], None, 0

    # once a read is mapped, record all mutations that were found in that read
    def find_snps(self, genome_window, mapped_idx, shift):
        min_len = min(len(self.sequence), len(genome_window))
        for i in range(min_len): 
            # check for SNPs by performing a letter-wise comparition b/w 
            # the read and the genome at the read's index
            if self.sequence[i] != genome_window[i]:
                mutation = f'>S{i+mapped_idx+shift} {genome_window[i]} {self.sequence[i]}'
                self.mutations.append(mutation)