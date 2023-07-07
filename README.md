# CS CM122 Project 1A: Read Mapper

For this project, I created a read mapper that mapped read to a reference genome using a Burrows Wheeler Transform. The output is list of all mutation that were contained within the reads, exluding sequencing errors.

The code for this project was run in the following format
    python3 project1a_main.py [reference genome] [reads]

For the results on Codalab, I ran
    python3 project1a_main.py project1a_10000_reference_genome.fasta project1a_10000_with_error_paired_reads.fasta