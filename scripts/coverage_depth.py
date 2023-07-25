#!/usr/bin/python 


import argparse
import os
import sys
import gzip
import re
from collections import defaultdict

parser = argparse.ArgumentParser(description = "\
DESCRIPTION: Script will take as input the files *.scrubbed_kmers.gz \n\
and *.kmer_hits.gz as well as the truth tabl (e.g., truth_vector_hits_only.txt \n\
the metagenomic targets (e.g., targets.txt) and the number of reads in each \n\
metagenome (e.g,. meta/metagenome_num_reads)\n\
It outputs the kmer_coverage, kmer_depth for each strain/metagenome combination\n\
")

# the input files
#parser.add_argument("--scrubbed_kmers_file", "-s", help = "File that contains the \
#set of informative kmers for a given strain based on scrubbing out frequent kmers in metagenomics and a large set of bacterial genomes", required = True)

#parser.add_argument("--metagenome_targets", "-e", help = "File that contains the \
#list of metagenomics files that were run against", required = True)

parser.add_argument("--kmer_hits_file", "-k", help = "File that contains the \
kmers with hits in each metagenomics file", required = True)

#parser.add_argument("--read_hits_file", "-r", help = "File that contains the \
#reads with hits in each metagenomics file", required = True)

parser.add_argument("--min_kmer_hits", "-m", help = "The minimum number of  \
kmers a read has matched to include the kmer from a read as a hit in the coverage. default 1; range (>=1)", required = False, default = 1, type=int)

#parser.add_argument("--truth_table_file", "-t", help = "Tab-delimited file with strain metagenome and 1 for truth  \
# (optional)", required = False)

parser.add_argument("--background_metagenomes_file", "-b", help = "file with background metagenome names  \
 (optional)", required = False)

#parser.add_argument("--num_metagenome_read_file", "-n", help = "Tab-delimited file metagenome file (PE1), number of reads, rounded grouping variable name for num reads (i.e., 10M, 0.1M, etc...) \
# (optional)", required = False)


#parser.add_argument("--independent", "-i", help = "Independently  \
#scrub metagenome and pangenome (without this flag, it scrubs them together so there is definitely -m min_fraction left)", action='store_true')

args = parser.parse_args()


# TO FIX: need to exclude commented out lines (startswith())
def count_total_kmers(scrubbed_kmer_file):
    line_count = 0
    with gzip.open(scrubbed_kmer_file, 'rt') as reader:
        for line in reader:
#            if not line.startswith('#'):
                line_count += 1
    
    #print(line_count)
    return line_count

def count_passed_kmers(kmer_hits_file, min_kmer_hits, kmer_read_hit_count):
    kmer_depth_count_by_metagenome = defaultdict(int)
    kmer_coverage_count_by_metagenome = defaultdict(int)
    unique_kmers_by_metagenome = defaultdict(int)
    kmer_total_evaluated_by_metagenome = defaultdict(int)
    read_total_evaluated_by_metagenome = defaultdict(int)
    genome_total_kmer = defaultdict(int)
    genome_total_informative_kmer = defaultdict(int)
    #all_unique_kmers = defaultdict(dict)

    passed_kmers = 0
    with gzip.open(kmer_hits_file, 'rt') as reader:
        for line in reader:
            if(not line.startswith('#')): # use # at the beginning to have comments and things
                content = line.rstrip('\n').split('\t')
                metagenomics_sample = os.path.basename(content[0])
                total_kmer_pe1 = int(content[1])
                total_informative_kmer_pe1 = int(content[2])
                total_kmer_pe2 = int(content[3])
                total_informative_kmer_pe2 = int(content[4])
                kmer_seq = content[5]

                total_informative_kmer = total_informative_kmer_pe1 + total_informative_kmer_pe2
                total_kmer = total_kmer_pe1 + total_kmer_pe2

                # TO FIX: this should be >= as it throws away kmers that occur in reads that have only 1 informative kmer (not the kmer's fault)
                # CORRECTION: actually I'm wrong as this is just operating on total_kmers
                if total_kmer > min_kmer_hits:
                     unique_str = metagenomics_sample + kmer_seq
#                     print(unique_str)
                     if unique_str not in unique_kmers_by_metagenome:
                         kmer_coverage_count_by_metagenome[metagenomics_sample] += 1
                         unique_kmers_by_metagenome[unique_str] = 1

                     kmer_depth_count_by_metagenome[metagenomics_sample] += 1
                   # if kmer_read_hit_count: # want to calculate the number of kmers hit on the read to get an average depth of the kmers across all reads
                    #    kmer_coverage_by_metagenome[metagenomics_sample] += total_informative_kmer # if a read had 10 informative kmers, then each kmer gets a "hit"
                    #else: # just want to count how many kmers were covered with passing reads
                     #   kmer_coverage_by_metagenome[metagenomics_sample] += 1
            else:
#                print(line)
                pieces = line.rstrip().split("\t")
                metagenomics_sample = os.path.basename(pieces[0])
                metagenomics_sample = re.sub('^#', "", metagenomics_sample)
                variable = pieces[1]
                value = int(pieces[2])
                if variable == "total_kmer_evaluated":
                    kmer_total_evaluated_by_metagenome[metagenomics_sample] = value
                elif variable == "total_reads_evaluated":
                    read_total_evaluated_by_metagenome[metagenomics_sample] = value
                elif variable == "total_genome_kmers":
                    genome_total_kmer[metagenomics_sample] = value
                elif variable == "total_genome_informative_kmers":
                    genome_total_informative_kmer[metagenomics_sample] = value

#                print("adding ", metagenomics_sample, " ", variable, " ", value)
                    

    # if we see the global statistics for a metagenome but not any informative kmer counts, it means the metagenome had no informative kmer counts, set to 0
    for metagenome in kmer_total_evaluated_by_metagenome:
        if not kmer_depth_count_by_metagenome[metagenome]:
            kmer_coverage_count_by_metagenome[metagenome] = 0
            kmer_depth_count_by_metagenome[metagenome] = 0
    
 
    
    #print(kmer_coverage_by_metagenome)
    return kmer_depth_count_by_metagenome, kmer_coverage_count_by_metagenome, kmer_total_evaluated_by_metagenome, read_total_evaluated_by_metagenome, genome_total_kmer, genome_total_informative_kmer

def get_background_meta(background_file):

    background_meta = {}
    with open(background_file, "rt") as reader:
        for line in reader:
            key = line.rstrip('\n')
            background_meta[key]=1
#            print("hit")

    return background_meta


def get_truth_table(truth_file):

    truth_table = {}
    with open(truth_file, "rt") as reader:
        for line in reader:
            content = line.rstrip('\n').split('\t')
            key = content[0] + "__" + content[1]
        #    print("key1 ", content[0], content[1], content[2], content[3], key)
            if int(content[2]) == 1:
                truth_table[key]="trueGenome"
            elif int(content[2]) == 2:
                truth_table[key]="trueMeta"
            elif int(content[2]) == 3:
                truth_table[key]="trueSpike"
            elif int(content[2]) == 10:
#                print('in type 10')
                truth_table[key]=content[3]


    return truth_table

def get_metagenomic_read_counts(metagenomic_read_count_file):

    count_table = {}
    with open(metagenomic_read_count_file, "rt") as reader:
        for line in reader:
            content = line.rstrip('\n').split('\t')
#            print("working on ", content[0], content[1], content[2])
           # count_table[content[0]] = [content[1], content[2]]
            count_table[content[0]] = [content[1], content[1]]
#                print("hit")
    return count_table
            
def base_metagenome_name(name):
#    print('name1', name)
    name = re.sub('_\d+M_PE\d+.fasta.gz$', "", name)
    name = re.sub('_\d+M_PE\d+.fastq.gz$', "", name)
    name = re.sub('_\d+K_PE\d+.fasta.gz$', "", name)
    name = re.sub('_\d+K_PE\d+.fastq.gz$', "", name)
    name = re.sub('_PE\d+.fasta.gz$', "", name)
    name = re.sub('_PE\d+.fastq.gz$', "", name)
    name = re.sub('_R\d+.fasta.gz$', "", name)
    name = re.sub('_R\d+.fastq.gz$', "", name)
    name = re.sub('.fasta.gz$', "", name)
    name = re.sub('.fastq.gz$', "", name)
#    print('name2', name)

    return name

def read_metagenome_files(metagenome_file_list):
    metagenome_files = {}
    with open(metagenome_file_list, "rt") as reader:
        for line in reader:
            clean = line.rstrip('\n')
            metagenome_files[clean] = 1
    return metagenome_files 

def main():
#    num_total_kmers = count_total_kmers(args.scrubbed_kmers_file)
    kmers_depth_count_in_metagenome, kmers_coverage_count_in_metagenome, kmer_total_evaluated_by_metagenome, read_total_evaluated_by_metagenome, genome_total_kmer, genome_total_informative_kmer = count_passed_kmers(args.kmer_hits_file, args.min_kmer_hits, True)

    strain_name = os.path.basename(args.kmer_hits_file)
    strain_name = re.sub(".kmer_hits.gz$", "", strain_name)
    species_pieces = strain_name.split('_')

    if (len(species_pieces)>1):
        species_name = species_pieces[0] + "_" + species_pieces[1]
    else:
        species_name = strain_name

    genus_name = species_pieces[0]


    # TO DO: parse the truth table and add that to the print out below (i.e., if it is a know true hit or not
#    truth_table = {}
    background_meta = {}
    if (args.background_metagenomes_file):
        background_meta = get_background_meta(args.background_metagenomes_file)


   # print("strain_name\tspecies_name\tgenus_name\tnum_informative_kmers\tmetagenome\tunique_observed_informative_kmers\ttotal_observed_informative_kmers\tkmer_coverage\tkmer_depth\tbase_metagenome\ttruth\tmetagenomic_reads\tmetagenome_reads_group\tkmer_depth_per_20B_kmer")
#    print("strain_name\tspecies_name\tgenus_name\tgenome_num_total_kmers\tgenome_num_informative_kmers\tmetagenome\tnum_metagenomic_reads\tnum_metagenome_kmers\tunique_observed_informative_kmers\ttotal_observed_informative_kmers\tkmer_coverage\tkmer_depth\tkmer_depth_per_20B_kmer\tbase_metagenome\tbackground")
    print("strain_name\tspecies_name\tgenus_name\tgenome_num_total_kmers\tgenome_num_informative_kmers\tmetagenome\tnum_metagenomic_reads\tnum_metagenome_kmers\tunique_observed_informative_kmers\ttotal_observed_informative_kmers\tkmer_coverage\tkmer_depth\tkmer_depth_per_20B_kmer\tbackground")
#    for metagenome in metagenome_targets:
    for metagenome in kmers_depth_count_in_metagenome:
        read_count = "NA"
        read_count_group = "NA"

#	print(metagenome)
        
        num_observed_kmers_in_metagenome = -1
        num_observed_unique_kmers_in_metagenome = -1
        num_evaluated_kmers_in_metagenome = -1
        num_reads_in_metagenome = -1
        num_genome_total_kmer = -1
        num_genome_total_informative_kmer = -1

#        print("on metagenome " , metagenome)
        if (metagenome in kmers_depth_count_in_metagenome):
            num_observed_kmers_in_metagenome = kmers_depth_count_in_metagenome[metagenome]
        if (metagenome in kmers_coverage_count_in_metagenome):
            num_observed_unique_kmers_in_metagenome = kmers_coverage_count_in_metagenome[metagenome]
        if (metagenome in kmer_total_evaluated_by_metagenome):
            num_evaluated_kmers_in_metagenome = kmer_total_evaluated_by_metagenome[metagenome]
        if (metagenome in kmer_total_evaluated_by_metagenome):
            num_reads_in_metagenome = read_total_evaluated_by_metagenome[metagenome]
        if (metagenome in genome_total_kmer):
            num_genome_total_kmer = genome_total_kmer[metagenome]
        if (metagenome in genome_total_informative_kmer):
            num_genome_total_informative_kmer = genome_total_informative_kmer[metagenome]

        kmer_coverage = num_observed_unique_kmers_in_metagenome / float(num_genome_total_informative_kmer)
        kmer_depth = num_observed_kmers_in_metagenome / float(num_genome_total_informative_kmer)
#	print('HERE c' + str(kmer_coverage) + ' d '+ str(kmer_depth) ,  num_observed_unique_kmers, num_observed_total_kmers, num_total_kmers)

        kmer_scale_constant = 2000000000
        if (num_evaluated_kmers_in_metagenome == 0):
            kmer_depth_scale = 0
        else:
            kmer_depth_default_kmer_scale_factor = kmer_scale_constant / float(num_evaluated_kmers_in_metagenome)
            kmer_depth_scale = kmer_depth * kmer_depth_default_kmer_scale_factor
#            print("scaling terms ", kmer_depth, kmer_scale_constant, num_evaluated_kmers, kmer_depth_default_kmer_scale_factor, kmer_depth_scale)


        background_status = 0
        if (metagenome in background_meta):
            background_status = 1

        print(strain_name + "\t" + species_name + "\t" + genus_name + "\t" + str(num_genome_total_kmer) + "\t" + str(num_genome_total_informative_kmer) + "\t" + metagenome + "\t" + str(num_reads_in_metagenome) + "\t" + str(num_evaluated_kmers_in_metagenome) + "\t" + str(num_observed_unique_kmers_in_metagenome) + "\t" + str(num_observed_kmers_in_metagenome) + "\t" + str(kmer_coverage) + "\t" + str(kmer_depth) + "\t" + str(kmer_depth_scale) + "\t" + str(background_status))

            

if __name__ == "__main__":
        main()
