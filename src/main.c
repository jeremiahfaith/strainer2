#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include "kseq.h"
#include "BIO_sequence.h"
#include "BIO_hash.h"
#include "genome_compare.h"
KSEQ_INIT(gzFile, gzread)

#define DEFAULT_SUBSAMPLE_SIZE 50000
#define DEFAULT_SEED 20
#define ALL_SEEDS 0
#define CLONE_MODE_SEEDS 50000
#define STRAIN_MODE_SEEDS 100000
#define STRAIN_MODE_THRESHOLD 0.05
#define CLONE_MODE_THRESHOLD 0.1
#define HYBRID_DEFAULT_THRESHOLD 0.1

/* 
	TO DO

	1) switch to tommy hash (fixed size) to see if faster

*/

void usage();

int main(int argc, char *argv[])
{
	BIO_hash seqHash;
	BIO_sequences seqs;
	int genome_length = 0;
	char *a_file = NULL;
	char *b_file = NULL; // standard file
	char *B_file = NULL; // file with a list of file names to compare with [so only one hash necessary]
	int seed = DEFAULT_SEED;
	int rapid_mode = ALL_SEEDS;
	int c;
	int print_header = 0;
	double threshold_for_fullmap = HYBRID_DEFAULT_THRESHOLD;
	int cloneMode = 0;
	int strainMode = 0;


        while ((c = getopt(argc, argv, "a:b:B:s:r:t:CSHhu")) != EOF)
                switch (c)
                {
                        case 'a': a_file    = strdup(optarg); break;
                        case 'b': b_file    = strdup(optarg); break;
                        case 'B': B_file    = strdup(optarg); break;
///                        case 's': seed      = strtol(optarg, NULL, 0); break;
                        case 's': seed      = atoi(optarg); break;
//                        case 'r': rapid_mode = DEFAULT_SUBSAMPLE_SIZE; break;
                        case 'r': rapid_mode = atoi(optarg); break;
                        case 't': threshold_for_fullmap = atof(optarg); break;
                        case 'H': print_header = 1; break;
                        case 'C': cloneMode = 1; break;
                        case 'S': strainMode = 1; break;
                        case 'u': usage(); break;
                        case 'h': usage(); break;
                        default: usage(); break;
                }




	if (!a_file || (!b_file && !B_file)) {
		usage();
		return 1;
	}
	
	if (print_header)
		printf("a_file\tb_file\thits\tmisses\tfrac\n");

	if (cloneMode) {
		rapid_mode = CLONE_MODE_SEEDS;
		threshold_for_fullmap = CLONE_MODE_THRESHOLD;
//		printf("clone mode %d\t%f\n", rapid_mode,threshold_for_fullmap);
	}
	if (strainMode) {
		rapid_mode = STRAIN_MODE_SEEDS;
		threshold_for_fullmap = STRAIN_MODE_THRESHOLD;
//		printf("strain mode %d\t%f\n", rapid_mode,threshold_for_fullmap);
	}
	if (cloneMode && strainMode) {
		fprintf(stderr, "Cannot run in clone mode and strain mode at same time (they are mutually exclusive)\n");
		usage();
		return 1;
	}

//		return 1;
//	}

	seqs = GEN_read_seq_file(a_file, &genome_length);
	//printf("genome length is %d\n", genome_length);
	seqHash = GEN_hash_sequences(seqs, seed, genome_length);
	//printf("hash done\n");

	if (b_file) { // only one vs one
		GEN_calculate_coverage(a_file, b_file, seed, seqHash, rapid_mode, threshold_for_fullmap);
	}
	else if (B_file) { // one vs many (can keep same hash)
		GEN_all_coverage(a_file, B_file, seed, seqHash, rapid_mode, threshold_for_fullmap);
	}
	
	

	// clean up
	BIO_destroySequences(seqs);
	BIO_destroyHash(seqHash);
	free(a_file);
	free(b_file);
	free(B_file);
	return 0;
}

void usage() {
	fprintf(stderr, "\n\nUsage: genome_compare -a <in1.fasta> -b <in2.fasta>\n");
	fprintf(stderr, "(alternative)         -a <in1.fasta> -B <file with multiple filenames>\n");
	fprintf(stderr, "                     [-s\tseed size (default=%d)]\n", DEFAULT_SEED);
	fprintf(stderr, "                     [-r\trapid mode (hashes entire reference compares with the first %d kmers in the query)]\n", DEFAULT_SUBSAMPLE_SIZE);
	fprintf(stderr, "                     [-t\tthreshold for fullmap; range (0.0 - 1.0); default %f (in rapid mode scores above this\n", HYBRID_DEFAULT_THRESHOLD);
	fprintf(stderr, "                        \tvalue are mapped again with the entire genome to get a more precise genome coverage score)]\n");
	fprintf(stderr, "\n\n\n");
	fprintf(stderr, "\t\tNote that rapid mode and the threshold are far more efficient when run with -B, as\n"); 
	fprintf(stderr, "\t\tthe large genome is hashed at the onset and the subsequent comparisions are extremely\n"); 
	fprintf(stderr, "\t\tfast if only using a subset of the genome\n");
}



