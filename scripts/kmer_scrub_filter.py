#!/usr/bin/python 

import argparse
import sys
import gzip
import re

parser = argparse.ArgumentParser(description = "\
DESCRIPTION: Script will take as input the output from the program kmer_scrub_count \n\
that lists every kmer in a strain of interest and the frequency of that kmer in a \n\
set of pangenomes and a set of metagenomes.\n\
With this program we then filter the set of kmers for the strain to select those that are most informative\n\
")

# the input file
parser.add_argument("--scrub_count_file", "-s", help = "This is the input file \
with the kmer counts to the pangenome and metagenomes", required = True)

parser.add_argument("--min_fraction", "-m", help = "The minimum fraction of  \
kmers to keep (so don't scrub everything); default 0.04; range (0.0-1.0)", required = False, default = 0.04, type=float)

parser.add_argument("--independent", "-i", help = "Independently  \
scrub metagenome and pangenome (without this flag, it scrubs them together so there is definitely -m min_fraction left)", action='store_true')

args = parser.parse_args()

def scrub_max_kmers(min_frac_to_keep, kmer_hash, total_kmers):
    min_kmer_count = -1
    fraction_kept = -1.0
    kmers_to_scrub = {} 
    total_kmers = float(total_kmers)

   # should add a quick version that first checks if the keys are fewer than threshold in which case it is waste of time to do the while loop
   # you would just return the keys

    # here we keep making our threshold more stringent until the number of kmers left drops below the threshold
    while(fraction_kept < min_frac_to_keep): 
        min_kmer_count+=1
        kmer_hits = 0

        for key,val in kmer_hash.items():
#            print(key + " " + str(val))
            if (val > min_kmer_count): 
                kmer_hits += 1

        fraction_kept = 1-(kmer_hits/total_kmers)
        sys.stderr.write("kept " + str(fraction_kept) + " with threshold " + str(min_kmer_count) + "\n")


    for key,val in kmer_hash.items():
        if (val > min_kmer_count): 
            kmers_to_scrub[key]=val

    sys.stderr.write("threshold was " + str(min_kmer_count) + " left with " + str(len(kmers_to_scrub)) + " out of " + str(total_kmers) + " that will be scrubbed\n")
    return kmers_to_scrub



def drug_scrub(drug_genome_hash, strain_hash, all_kmers):
    strain_hash_scrubbed = strain_hash
    for key in drug_genome_hash:
        if key in strain_hash_scrubbed:
            del strain_hash_scrubbed[key]

    return strain_hash_scrubbed



def independent_scrub(min_fraction, pangenome_hash, metagenome_hash, strain_hash, all_kmers):
    strain_hash_scrubbed = strain_hash
    kmers_to_scrub_pan = scrub_max_kmers(min_fraction, pangenome_hash, all_kmers)
    kmers_to_scrub_meta = scrub_max_kmers(min_fraction, metagenome_hash, all_kmers)

    for key in kmers_to_scrub_pan:
        if key in strain_hash_scrubbed:
            del strain_hash_scrubbed[key]
    for key in kmers_to_scrub_meta:
        if key in strain_hash_scrubbed:
            del strain_hash_scrubbed[key]
    
    return strain_hash_scrubbed


def joint_scrub(min_fraction, pangenome_hash, metagenome_hash, strain_hash, all_kmers, num_drug_scrubbed): 

    # convert from counts to percentages
    metagenome_sum = 0
    for key in metagenome_hash:
        metagenome_sum += metagenome_hash[key]

    for key in metagenome_hash:
#        print("mN" + str(metagenome_hash[key]))
        metagenome_hash[key]/=float(metagenome_sum)
#        print("mP" + str(metagenome_hash[key]))
#    print("meta " + str(metagenome_sum))

    pangenome_sum = 0
    for key in pangenome_hash:
        pangenome_sum += pangenome_hash[key]

    for key in pangenome_hash:
#        print("pN " + str(pangenome_hash[key]))
        pangenome_hash[key]/=float(pangenome_sum)


    # combine the two hashes so we remove the kmers in order of frequency
    for key in strain_hash:
        strain_hash[key]=0
        if key in metagenome_hash and metagenome_hash[key] > strain_hash[key]:
            strain_hash[key] = metagenome_hash[key]
        if key in pangenome_hash and pangenome_hash[key] > strain_hash[key]:
            strain_hash[key] = pangenome_hash[key]

#        print("pP " + str(pangenome_hash[key]))
#    print("pan " + str(pangenome_sum))
    pairs = sorted(strain_hash.items(), key=lambda x: x[1], reverse=True)


    # remove them from biggest to smallest (most frequent to least frequent)
    num_scrubbed = float(num_drug_scrubbed)
    scrub = {}
    for pair in pairs:
        if (1-((num_scrubbed+1)/all_kmers)) > min_fraction:
#            print(pair)
            scrub[pair[0]] = pair[1]
            num_scrubbed += 1.0

            if pair[0] in strain_hash:
                del strain_hash[pair[0]]
#                print("removing " + key)
#        else:
#            sys.stderr.write("scrubbed enough " + str(num_scrubbed/all_kmers) + "\n")

    return strain_hash
   


def main():
    if args.min_fraction < 0.0 or args.min_fraction > 1.0:
        sys.stderr.write("error --min_fraction (-m) must be between 0.0 and 1.0 (" + args.min_fraction +  ")\n")

    strain_hash = {}
    metagenome_hash = {}
    pangenome_hash = {}
    drug_genome_hash = {}
    drug_filter = 0
    all_kmers = 0
    with gzip.open(args.scrub_count_file, 'rt') as reader:
        for line in reader:
            if(not line.startswith('#')): # use # at the beginning to have comments and things
                content = line.rstrip('\n').split('\t')
                key = content[0]
                all_kmers += 1
                content[1] = int(content[1])
                content[2] = int(content[2])
                content[3] = int(content[3])
                strain_hash[key] = content[1]
                

                if (content[2] > 0):
                    pangenome_hash[key] = content[2]
                if (content[3] > 0):
                    metagenome_hash[key] = content[3]

                if (len(content) == 5):
                    drug_filter = 1
                    content[4] = int(content[4])        
                    if (content[4] > 0):
                        drug_genome_hash[key] = content[4]
                #strain_hash[key] = content[1]
                #print("A:" + content[0] + " B:" + content[1] + " C:" + content[2])

    print("#total kmers is strain:" + str(all_kmers) + "," + str(len(strain_hash)) + " pangenome: " + str(len(pangenome_hash)) + " metagenome: " + str(len(metagenome_hash)))


    drug_scrubbed = 0 # how many kmers are scrubbed because they overlap between strains in the same drug
    if (drug_filter): 
        print("#total kmers cross drug:" + str(len(drug_genome_hash)))
        strain_hash = drug_scrub(drug_genome_hash, strain_hash, all_kmers)
        fraction_kmer_remaining_post_drug_scrub = float(len(strain_hash)/float(all_kmers));
        drug_scrubbed = all_kmers - len(strain_hash)
        print("#fraction kmers remaining drug post scrub:" + str(fraction_kmer_remaining_post_drug_scrub))
        print("#drug_scrubbed kmers:" + str(drug_scrubbed))

        # don't go on if not enough kmers left after the drug scrub
        if (fraction_kmer_remaining_post_drug_scrub < args.min_fraction * 2):
            raise Exception('ERROR: too few kmers remain after drug scrub. Are your drug strains too similar?')

    

    if args.independent:
        strain_hash_scrubbed = independent_scrub(args.min_fraction, pangenome_hash, metagenome_hash, strain_hash, all_kmers)
    else:
        strain_hash_scrubbed = joint_scrub(args.min_fraction, pangenome_hash, metagenome_hash, strain_hash, all_kmers, drug_scrubbed) 


    print("#post scrub kmers " + str(len(strain_hash_scrubbed)) + " out of " + str(all_kmers))

    for key in strain_hash_scrubbed:
        print(key)
            

if __name__ == "__main__":
        main()
