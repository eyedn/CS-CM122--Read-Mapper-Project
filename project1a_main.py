###############################################################################
#   Aydin Karatas
#   CS 122 S23
#   Project 1a
#   project1a_main.py
###############################################################################
import project1a_classes as cl
import project1a_functions as fx
from sys import argv

def main():
    genome_file = argv[1]
    reads_file = argv[2]
    # create genome objcet and dictionary of reads objects
    my_genome = cl.Genome(genome_file)
    my_reads = fx.create_reads_dict(reads_file, 10, 7)
    # map read objects to genome object
    fx.map_reads(my_genome, my_reads, 2)
    # consensus algorithms:
    # only consider mutations that appear more than once to ignore read errors
    my_mutations = fx.concensus_mutations(my_reads, 2)
    for mut in my_mutations:
        print(mut)


if __name__ == "__main__":
    main()