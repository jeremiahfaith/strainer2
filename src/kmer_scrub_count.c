#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include "kseq.h"
#include "BIO_sequence.h"
#include "BIO_hash.h"
#include "genome_compare.h"
#include <inttypes.h>
KSEQ_INIT(gzFile, gzread)


/* 
	TO DO

	1) switch to tommy hash (fixed size) to see if faster


	TIMINGS
	took 16min on 1255 non redundant genomes in latest version of the culture database (Sept 2020)
	took 68min on 1255 non redundant genomes in latest version of the culture database (Sept 2020) and FOCUS metagenomes

*/

void usage();
void calculate_coverage(const char *a_file, const char *file, const int seed, BIO_hash h);
void all_coverage(const char *a_file, const char *B_file, const int seed, BIO_hash seqHash);
void print_hash_counts(BIO_hash seqHash, const char *C_file);

int main(int argc, char *argv[])
{
	BIO_hash seqHash;
//	BIO_sequences seqs;
//	int genome_length = 0;
	char *A_file = NULL; // file with a list of genome file names to count scrubbing from
	char *B_file = NULL; // file with a list of metagenome file names to count scrubbing from
	char *C_file = NULL; // file with a list of genome file names of the drug genomes
	char *r_file = NULL; // reference genome file
	char *p_file = NULL; // optional progress file
	const int seed = 31; // using the same number that Varun used for consistency
	const int default_hash_val = 1;
	const int default_hash_increment = 1;
	//const int size_of_hash_vec = 3; // place for original count, pangenome count, metagenome count
	const int size_of_hash_vec = 4; // place for original count, pangenome count, metagenome count, drug count
	int hash_index = 0; // to switch from original count, pangenome count, metagenome count
//	const int seed = 64;
//	const int seed = 96;
	int c;
	FILE *fout;
	FILE *progress = NULL;
	int write_dist = 0;

        while ((c = getopt(argc, argv, "A:B:C:r:p:Hhud")) != EOF)
                switch (c)
                {
              //          case 'a': a_file    = strdup(optarg); break;
                        case 'A': A_file    = strdup(optarg); break;
                        case 'B': B_file    = strdup(optarg); break;
                        case 'C': C_file    = strdup(optarg); break;
                        case 'r': r_file    = strdup(optarg); break;
                        case 'p': p_file    = strdup(optarg); break;
               //         case 'b': b_file    = strdup(optarg); break;
                //        case 'B': B_file    = strdup(optarg); break;
                        case 'd': write_dist = 1; break;
                        case 'u': usage(); break;
                        case 'h': usage(); break;
                        default: usage(); break;
                }



//	if ((!a_file && !A_file) || (!b_file && !B_file)) {
	if (!r_file || !A_file || !B_file) {
		usage();
		return 1;
	}


	if (p_file != NULL) {
		progress = fopen(p_file, "w");
		if (progress == NULL) {
			fprintf(stderr, "could not open progress file %s\n", p_file);
			exit(EXIT_FAILURE);
		}
		fprintf(progress, "adding kmer counts for:\n");
	}

        seqHash = BIO_initHash(DEFAULT_GENOME_HASH_SIZE);

	GEN_hash_sequences_set_count_vec(r_file, seed, seqHash, default_hash_val, default_hash_increment, 0, size_of_hash_vec);
	GEN_all_kmer_counts(A_file, seed, seqHash, 1, progress); // store this pangenome result in the second column 
	GEN_all_kmer_counts(B_file, seed, seqHash, 2, progress); // store this metagenome result in the third column 

	if (C_file) 
		GEN_all_kmer_counts_skip_file(C_file, r_file, seed, seqHash, 3, progress); // skip the hashing of the reference file



	print_hash_counts(seqHash, C_file);

	

//	if (b_file) // only one vs one
//	printf("checked %" PRIu64 " seeds\n", metagenome_seed_count);
/*
	else if (B_file) // one vs many (can keep same hash)
		all_coverage(a_file, B_file, seed, seqHash);
*/
	
	

	// clean up
	BIO_destroyHashD(seqHash);
//	free(a_file);
//	free(b_file);
	free(A_file);
	free(B_file);
	free(C_file);
	free(r_file);
	free(p_file);
	if (progress != NULL)
		fclose(progress);
//	free(B_file);
	return 0;
}

void usage() {
	fprintf(stderr, "Usage: kmer_scrub_count -r <reference genome>  -A <file with multiple genome filenames> -B <file with multiple metagenome filenames> -C <(optional) file with multiple genome filenames of drug strains> -p [progress output file, optional]\n");
//	fprintf(stderr, "                         -r [reference genome; if not provided writes pangenome for all genomes in -A]\n");
//	fprintf(stderr, "                         -d [print all hash counts (e.g., to plot distribution across entire pangenome]\n");
//	fprintf(stderr, "(alternative)         -a <in1.fasta> -B <file with multiple filenames>\n");
}


void print_hash_counts(BIO_hash seqHash, const char *C_file) {
	char **allKeys = BIO_getHashKeys(seqHash);
	int hash_size = BIO_getHashSize(seqHash);
	unsigned int i;
        unsigned int *counts = NULL;

//	printf("have %d keys\n", hash_size);
	//printf("#kmer\treference_count\tpangenome_count\tmetagenome_count\n");
	printf("#kmer\treference_count\tpangenome_count\tmetagenome_count\tdrug_count\n");
	for (i=0; i<hash_size; i++) {
		counts = (unsigned int*)BIO_searchHash(seqHash,allKeys[i]);
		
		if (C_file) {
			printf("%s\t%d\t%d\t%d\t%d\n", allKeys[i], counts[0], counts[1], counts[2], counts[3]);
		}
		else {
			printf("%s\t%d\t%d\t%d\n", allKeys[i], counts[0], counts[1], counts[2]);
		}
	}


	BIO_destroyHashKeys(allKeys);
}
