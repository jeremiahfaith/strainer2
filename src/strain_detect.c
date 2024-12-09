#include "zlib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "kseq.h"
#include "BIO_sequence.h"
#include "BIO_hash.h"
#include "genome_compare.h"
#include <inttypes.h>
#define KMER_TYPE 0
#define INFORMATIVE_KMER_MATCH_PE1 1
#define TOTAL_KMER_MATCH_PE1 2
#define INFORMATIVE_KMER_MATCH_PE2 3
#define TOTAL_KMER_MATCH_PE2 4
#define BACKGROUND_KMER_COUNT 5
#define INFORMATIVE_KMER 2
#define NON_INFORMATIVE_KMER 1
#define NOT_PAIRED_END 0
#define IS_PAIRED_END 1
#define IS_PAIRED_END_INTERLEAVE 2
#define UNKNOWN_FILE_TYPE -1
//#define NO_GZIP_OUTPUT 1

// these could actually be stored as a char or short to save space because the counts is for a PE150nt read [thus < 256] but 
// would need to remember later if used on a machine with larger read
// first column is for the type (informative_kmer or non_informative_kmer); second column is to count max informative_kmer hits in a read
// and third column is to count the number of total kmer hits in a read (for species level information)
// will ideally make this size 5 if we switch to do for PE1 and PE2
#define HASH_VECTOR_SIZE 6
//#define HASH_VECTOR_SIZE 3
KSEQ_INIT(gzFile, gzread)


/* 
	TO DO

        1) switch to consider non-informative kmers (to look for high quality match to the species) and informative kmers (to look for 
        informative information about the strain). The hope is that this will clean up a little the precision because we might find 
        sufficient informative kmers where there is not a good species match; for now doing wih PE1 only but this should get substantially 
        better if considering PE2 and PE1 for the non-informative kmers as that will be more stringent for species hit

	2) if we don't need to store the top overall hit for a kmer anymore, we could remove this vector of numbers as it is no longer used
//		count[TOTAL_KMER_MATCH_PE1]=read_kmer_hits;
//		count[TOTAL_KMER_MATCH_PE2]=read_kmer_hits_pe2;
//		count[INFORMATIVE_KMER_MATCH_PE1]=read_informative_kmer_hits;
//		count[INFORMATIVE_KMER_MATCH_PE2]=read_informative_kmer_hits_pe2;


*/

void usage();
unsigned int hash_scrubbed_kmers(const char *filename, BIO_hash h, const int seed);
void quantify_hits_all_files(const char *file_of_filenames, const char *file_PE1, const char *file_PE2, BIO_hash h, const int seed, const char *outfile, int is_PE);
void quantify_hits_PE(const char *PE1_file, const char *PE2_file, BIO_hash h, const int seed, FILE *out, gzFile gzout, int hash_size, const int is_PE, const unsigned int genome_kmers, const unsigned int genome_informative_kmers);
void background_filter(const char *background_file, const double fraction_to_remove, const unsigned int num_inform_kmer, BIO_hash h, const int seed);
unsigned int get_file_type(const char *filetype);
int compare (const void *a, const void *b);
int kmer_removed(int max_kmer_to_keep, unsigned int *informative_kmer_background_counts, int vec_size);

int main(int argc, char *argv[])
{
	BIO_hash seqHash;
//	BIO_sequences seqs;
//	int genome_length = 0;
	char *a_file = NULL;
	char *r_file = NULL;
	char *b_file = NULL; // standard file
	char *b_file2 = NULL; // standard file
	char *A_file = NULL; // file with a list of file names to compare with [so only one hash necessary]
	char *B_file = NULL; // file with a list of file names to compare with [so only one hash necessary]
	char *seq_file_type = NULL; // PE, SE, IPE
	char *background_file = NULL; // file with a list of file names to metagenomic files that represent background
	char *kmer_outfile = NULL;
	int is_paired_end = NOT_PAIRED_END;
	unsigned int num_informative_kmers;
//	const int seed = 32;
	const int seed = 31;
	int c;
	FILE *fout;
	int initial_size = 1000000;
	double fraction_background_to_remove = 0.5; // default to removing 1/2 of the background signal

        while ((c = getopt(argc, argv, "g:r:a:A:b:c:B:S:M:o:t:Hhuspn")) != EOF)
                switch (c)
                {
                        case 'a': a_file    = strdup(optarg); break;
                        case 'A': A_file    = strdup(optarg); break;
                        case 'b': b_file    = strdup(optarg); break;
                        case 'c': b_file2    = strdup(optarg); break;
                        case 'B': B_file    = strdup(optarg); break;
                        case 'r': r_file    = strdup(optarg); break;
                        case 'g': background_file    = strdup(optarg); break;
                        case 'o': kmer_outfile = strdup(optarg); break;
                        case 'n': is_paired_end = NOT_PAIRED_END; break;
                        case 't': seq_file_type = strdup(optarg); break;
                        case 'u': usage(); break;
                        case 'h': usage(); break;
                        default: usage(); break;
                }


//	if ((!a_file && !A_file) || (!b_file && !B_file)) {
	if (!a_file || !kmer_outfile || !r_file) {
		usage();
		exit(1);
	}
	if (!b_file && !B_file) {
		usage();
		exit(1);
	}

	if (seq_file_type != NULL) {
		is_paired_end = get_file_type(seq_file_type);

		if (is_paired_end == UNKNOWN_FILE_TYPE) {
			printf("unknown filetype specification. allowed are SE, PE, PEI\n\n");
			usage();	
			exit(1);
		}
	}
	
	if (b_file && is_paired_end == IS_PAIRED_END) {
		if (!b_file2) {
			printf("commandline PE mapping requires two files (-b [file1] and -c [file2])\n\n");
			usage();
			exit(1);
		}
	}
	if (b_file != NULL && B_file != NULL) {
		printf("cannot have -B flag and -b flag\nEither have a file with metagenomics files to be detect the strain in or specify one metagenomic file to detect the strain in\n");
		usage();
		exit(1);
	}
	

	seqHash = BIO_initHash(DEFAULT_GENOME_HASH_SIZE);
//	printf("reading reference hash from genome sequence\n");
        GEN_hash_sequences_set_count_vec(r_file, seed, seqHash, NON_INFORMATIVE_KMER, 0, 0, HASH_VECTOR_SIZE);
	num_informative_kmers = hash_scrubbed_kmers(a_file, seqHash, seed);

	if (background_file) 
        	background_filter(background_file, fraction_background_to_remove, num_informative_kmers, seqHash, seed);


	quantify_hits_all_files(B_file, b_file, b_file2, seqHash, seed, kmer_outfile, is_paired_end);


	BIO_destroyHashD(seqHash);
	free(a_file);
	free(b_file);
	free(b_file2);
	free(A_file);
	free(B_file);
	free(background_file);
	free(seq_file_type);
	return 0;
}

void background_filter(const char *background_file, const double fraction_to_remove, const unsigned int num_inform_kmer, BIO_hash h, const int seed) {
	unsigned int kmer_to_keep = (int) (num_inform_kmer * fraction_to_remove);
	unsigned int *informative_kmer_background_counts = (unsigned int*)calloc(num_inform_kmer, sizeof(unsigned int));
	char **allKeys = BIO_getHashKeys(h);
	int hash_size = BIO_getHashSize(h);
	int kmer_count = 0;
        unsigned int *count = NULL;
	int i;
	int max_kmer_to_keep = 1;
	int actual_keep = 0;
	
	

	printf("#removing %f proportion of %s kmers; informative %d keep at least %d\n", fraction_to_remove, background_file, num_inform_kmer, kmer_to_keep);


//	printf("reading bkgnd kmer\n");	
	GEN_all_kmer_counts(background_file, seed, h, BACKGROUND_KMER_COUNT, NULL); // store this pangenome result in the second column 
//	printf("finished bkgd kmer read\n");
	for (i=0; i<hash_size; i++) {
		count = (unsigned int*)BIO_searchHash(h,allKeys[i]);
		if (count[KMER_TYPE] == INFORMATIVE_KMER) {
			if (kmer_count >= num_inform_kmer) {
				fprintf(stderr, "Error: too many background kmers\n");
				exit(1);
			}

			informative_kmer_background_counts[kmer_count] = count[BACKGROUND_KMER_COUNT];
//			if (informative_kmer_background_counts[kmer_count]) 
//				printf("back %d\n", informative_kmer_background_counts[kmer_count]);
			kmer_count++;
		}
		
	}

	qsort(informative_kmer_background_counts, num_inform_kmer, sizeof(unsigned int), compare);

	/// now go through and see the count at kmer_to_keep
	
//	for (i=0; i<num_inform_kmer; i++) {
//		printf("%d\t%d\n", i, informative_kmer_background_counts[i]);
///	}

	//printf("have %d informative\n", kmer_count);
	// now sort the background kmers to see what value we remove above; when there are alot of ties we have to be careful to not remove too many; guess will just remove in kmer order
	// must have at least count of 1 (the detault) so if above that we move the threshold up otherwise we remove all the kmers >= 1 
	
	if (informative_kmer_background_counts[kmer_to_keep-1] > max_kmer_to_keep) 
		max_kmer_to_keep = informative_kmer_background_counts[kmer_to_keep-1]; 

//	kmer_to_keep = 500; just for debugging making more stringent threshold to make sure works
	qsort(informative_kmer_background_counts, num_inform_kmer, sizeof(unsigned int), compare);
	while (kmer_removed(max_kmer_to_keep, informative_kmer_background_counts, num_inform_kmer) > kmer_to_keep) {
		max_kmer_to_keep++; // keep trying higher kmer background thresholds until you don't remove too many kmers
	}
	

	// now we want to run through and designate as uninformative all the kmers where the background is above the threshold
	kmer_count = 0;
	for (i=0; i<hash_size; i++) {
		count = (unsigned int*)BIO_searchHash(h,allKeys[i]);
		if (count[KMER_TYPE] == INFORMATIVE_KMER) {
			if (count[BACKGROUND_KMER_COUNT] >= max_kmer_to_keep) { // removing this kmer from the informative list
				count[KMER_TYPE] = NON_INFORMATIVE_KMER;
				kmer_count++;
			}
		}
		
	}
	printf("#final_threshold %d removes %d background kmers %d removed\n",max_kmer_to_keep, kmer_removed(max_kmer_to_keep, informative_kmer_background_counts, num_inform_kmer), kmer_count); 


	/// now go through and see the count at kmer_to_keep
	




	free(informative_kmer_background_counts);
        BIO_destroyHashKeys(allKeys);
}

int kmer_removed(int max_kmer_to_keep, unsigned int *informative_kmer_background_counts, int vec_size) {
	int i;
	int count = 0;

	for (i=0; i<vec_size; i++) {
		if (informative_kmer_background_counts[i] >= max_kmer_to_keep)
			count++;
	}

	return count;
}


int compare (const void *a, const void *b) {
//	unsigned int *x = (unsigned int *) a;
//	unsigned int *y = (unsigned int *) b;

//	return *x - *y;
	return (*(unsigned int *)b - *(unsigned int*)a);
}

void quantify_hits_all_files(const char *file_of_filenames, const char *file_PE1, const char *file_PE2, BIO_hash h, const int seed, const char *outfile, int is_PE) {
	FILE *fp;
	char *line1 = NULL;
	size_t len1 = 0;
	ssize_t read_s1;
	char *line2 = NULL;
	char *file1 = NULL;
	char *file2 = NULL;
	size_t len2 = 0;
	ssize_t read_s2;
	char *pos;
	char *token;
	char **allKeys = BIO_getHashKeys(h);
	int hash_size = BIO_getHashSize(h);
	unsigned int total_genome_kmers = hash_size;
	unsigned int total_genome_informative_kmers = 0;
	int i;
        unsigned int *count = NULL;
	gzFile gzout = NULL;
	FILE *out = NULL;


	for (i=0; i<hash_size; i++) {
		count = (unsigned int*)BIO_searchHash(h,allKeys[i]);
		if (count[KMER_TYPE] == INFORMATIVE_KMER) {
			total_genome_informative_kmers++;
		}
	}

	#ifdef NO_GZIP_OUTPUT
		out = fopen(outfile, "w");
		if (out == NULL) {
			fprintf(stderr, "could not open *out file outfile %s in quantify_hits_all_files()\n", outfile);
			exit(EXIT_FAILURE);
       		}
	#else
		gzout = gzopen(outfile, "wb9");
		if (gzout == NULL) {
			fprintf(stderr, "could not open *gzout file outfile %s in quantify_hits_all_files()\n", outfile);
			exit(EXIT_FAILURE);
       		}
	#endif

	if (file_of_filenames != NULL) {
		fp = fopen(file_of_filenames, "r");
		if (fp == NULL) {
			fprintf(stderr, "could not read file file_of_filenames %s in quantify_hits_all_files()\n", file_of_filenames);
			exit(EXIT_FAILURE);
		}

		while ((read_s1 = getline(&line1, &len1, fp)) != -1) {
			is_PE = -1;

			if ((pos=strchr(line1, '\n')) != NULL)
				*pos = '\0';

//			printf("here with read_s1\n");
	
			token = strtok(line1, "\t");
//			printf("file type is %d\n", is_PE);
			is_PE = get_file_type(token);
//			printf("file type is %d\n", is_PE);

			if (is_PE == UNKNOWN_FILE_TYPE) {
				printf("unknown file type skipping line (%s)\n", token);
			}
			else {
				file1 = strtok(NULL, "\t");
				if (file1 == NULL) {
					printf("ERROR: no first file specified for %s\n", line1);
				}
				else {
					if (is_PE == NOT_PAIRED_END) {
//						printf("running SE with %s\t%s\n", file1, line1);
						quantify_hits_PE(file1, NULL, h, seed, out, gzout, hash_size, is_PE, total_genome_kmers, total_genome_informative_kmers);
					}
					else if (is_PE == IS_PAIRED_END_INTERLEAVE) {
//						printf("running PEI with %s\t%s\n", file1, line1);
						quantify_hits_PE(file1, NULL, h, seed, out, gzout, hash_size, is_PE, total_genome_kmers, total_genome_informative_kmers);
					}
					else if (is_PE == IS_PAIRED_END) {
						file2 = strtok(NULL, "\t");
						if (file2 == NULL) {
							printf("ERROR: no second file specified for PE: %s\n", line1);
						}
						else {
//							printf("running PE with %s\t%s\t%s\n", file1, file2, line1);
							quantify_hits_PE(file1, file2, h, seed, out, gzout, hash_size, is_PE, total_genome_kmers, total_genome_informative_kmers);
						}
					}
				}
			}
/*
			if (is_PE) {
				read_s2 = getline(&line2, &len2, fp);
				if (read_s2 == -1)  {
					printf("ERROR: must be an even number of files one for each end of the paired-end read\n");
				}
				if ((pos=strchr(line2, '\n')) != NULL) {
					*pos = '\0';
				}
			}
*/

//			printf("%s\t%s\n", line1, line2);
//			quantify_hits_PE(line1, line2, h, seed, out, allKeys, hash_size, is_PE);
		}
		fclose(fp);
	}
	else {
		quantify_hits_PE(file_PE1, file_PE2, h, seed, out, gzout, hash_size, is_PE, total_genome_kmers, total_genome_informative_kmers);
	}

	#ifdef NO_GZIP_OUTPUT
	fclose(out);
	#else
	gzclose(gzout);
	#endif
	free(line1);
	free(line2);
        BIO_destroyHashKeys(allKeys);
}

/* this is where we spend most of the CPU time; have done quite a lot of optimization to speed up; major bottleneck is the hash; trying to improve that next... */
void quantify_hits_PE(const char *PE1_file, const char *PE2_file, BIO_hash h, const int seed, FILE *out, gzFile gzout, int hash_size, const int is_PE, const unsigned int genome_kmers, const unsigned int genome_informative_kmers) {
	gzFile fp;
	gzFile fp2 = NULL;
	kseq_t *seq, *seq2;
	int l, l2;
	unsigned int i;
//	char *seedStrRevComp = (char*)malloc(sizeof(char) * seed+1);
	char *seedStrRevComp = (char*)malloc(sizeof(char) * 10000);
	char *sequence_PE1_copy = (char*)malloc(sizeof(char) * 10000);
	int sequence_PE1_length = 0;
	char *orientStr;
	char temp_nuc;
	char *seed_seq;
        unsigned int *count = NULL;
	int read_kmer_hits = 0;
	int read_informative_kmer_hits = 0;
	int min_hits_for_informative_read = 1;
	int read_kmer_hits_pe2 = 0;
	int read_informative_kmer_hits_pe2 = 0;
	int min_hits_for_good_match = 1;
	int total_read_kmer_hits = 0;
	int total_read_informative_kmer_hits = 0;
	long long unsigned int total_kmers_evaluated = 0;
	long long unsigned int total_reads_evaluated = 0;
	int has_N;
//	int min_hits_for_informative_read = 5;
//	int min_hits_for_good_match = 5;
//	int min_hits_for_informative_read = 10;
//	int min_hits_for_good_match = 10;


	fp = gzopen(PE1_file, "r");
	if (fp == NULL) {
		fprintf(stderr, "could not read file (read1) %s in quantify_hits_PE() (error: %s)\n", PE1_file, strerror(errno));
		exit(EXIT_FAILURE);
	}
	seq = kseq_init(fp);


	if (is_PE) {
		if (is_PE == IS_PAIRED_END) {
			fp2 = gzopen(PE2_file, "r");
			if (fp2 == NULL) {
				fprintf(stderr, "could not read file (read2) is_PE %s in quantify_hits_PE() (error: %s)\n", PE2_file), strerror(errno);
				exit(EXIT_FAILURE);
			}
			seq2 = kseq_init(fp2);
		}
		else if (is_PE == IS_PAIRED_END_INTERLEAVE) { // for interleaved we just keep reading the same file
			seq2 = seq;
		}
	}

//	fprintf(stderr, "reading %s to check in hash %d\n", PE1, BIO_getHashSize(h));

	char *seed_seq_rc;
	while (l = kseq_read(seq) >= 0) {
		if (seq->seq.l >= seed) { // skip if a read has been trimmed shorter than the seed length
			total_reads_evaluated++;
			read_kmer_hits = 0;
			read_informative_kmer_hits = 0;
			sequence_PE1_length = seq->seq.l;
			BIO_stringToUpper(seq->seq.s); // keep same case

			seed_seq = seq->seq.s;

			// need to copy this stuff because in interleaved files it is forgotten when I read the next one
			strcpy(sequence_PE1_copy, seq->seq.s);

			// make the reverse complement now so that we can move along the string using pointers and avoid complementing substrings over and over
			strcpy(seedStrRevComp, seq->seq.s);
			BIO_reverseComplement(seedStrRevComp);
			seed_seq_rc = &seedStrRevComp[sequence_PE1_length-seed];

			/* do PE1 */

			has_N = contains_N(seed_seq);

			for (i = 0; i<seq->seq.l - seed+1; i++) { // for each possible seed position
				temp_nuc = seed_seq[seed];
				seed_seq[seed] = '\0';
				seed_seq_rc[seed] = '\0';
//				printf("%s\n%s\n\n", seed_seq, seedStrRevComp);
				//if (total_reads_evaluated == 10) { exit(0);}
                                if (strcmp(seed_seq, seed_seq_rc) > 0) //  orient the strings doing the upfront revcomp with pointer shifts is slightly faster than the orient_string function
                                        orientStr = seed_seq;
                                else
                                        orientStr = seed_seq_rc;
	
//				orientStr = orient_string(seed_seq, seedStrRevComp, seed);
				if (!has_N || !contains_N(orientStr)) {
					count = (unsigned int*)BIO_searchHash(h,orientStr);
					if (count != NULL) {
//						*count+=1; // increment the n-mer
						read_kmer_hits++; // all hits
						if (count[KMER_TYPE] == INFORMATIVE_KMER)
							read_informative_kmer_hits++; // informative hits
					}
				}
	
				seed_seq[seed] = temp_nuc;
				seed_seq++;
				seed_seq_rc--;
				total_kmers_evaluated++;
			}
		}

		/* do PE2 */
		if (is_PE) {
			l2 = kseq_read(seq2);
			if (seq2->seq.l >= seed) { // skip if a read has been trimmed shorter than the seed length

				read_kmer_hits_pe2 = 0;
				read_informative_kmer_hits_pe2 = 0;
				if (l2 < 0) {
					fprintf(stderr, "reached end of PE2 (%s) before end of PE1 (%s), check that file names are correct\n", PE2_file, PE1_file);
					exit(EXIT_FAILURE);
				}
				BIO_stringToUpper(seq2->seq.s); // keep same case
				seed_seq = seq2->seq.s;
				has_N = contains_N(seed_seq);

				// make the reverse complement now so that we can move along the string using pointers and avoid complementing substrings over and over
				strcpy(seedStrRevComp, seq2->seq.s);
				BIO_reverseComplement(seedStrRevComp);
				seed_seq_rc = &seedStrRevComp[seq2->seq.l - seed];

				for (i = 0; i<seq2->seq.l - seed+1; i++) { // for each possible seed position
					temp_nuc = seed_seq[seed];
					seed_seq[seed] = '\0';
					seed_seq_rc[seed] = '\0';
		
				//	orientStr = orient_string(seed_seq, seedStrRevComp, seed);
                                	if (strcmp(seed_seq, seed_seq_rc) > 0) //  orient the strings doing the upfront revcomp with pointer shifts is slightly faster than the orient_string function
                                	        orientStr = seed_seq;
                                	else
                                	        orientStr = seed_seq_rc;

					if (!has_N || !contains_N(orientStr)) {
						count = (unsigned int*)BIO_searchHash(h,orientStr);
						if (count != NULL) {
//							*count+=1; // increment the n-mer
							read_kmer_hits_pe2++; // all hits
							if (count[KMER_TYPE] == INFORMATIVE_KMER)
								read_informative_kmer_hits_pe2++; // informative hits
						}
					}
		
					seed_seq[seed] = temp_nuc;
					seed_seq++;
					seed_seq_rc--;
					total_kmers_evaluated++;
				}
			}
		}

		total_read_kmer_hits = read_kmer_hits + read_kmer_hits_pe2;
		total_read_informative_kmer_hits = read_informative_kmer_hits + read_informative_kmer_hits_pe2;

//		if ((read_kmer_hits >= min_hits_for_good_match || read_kmer_hits_pe2 >= min_hits_for_good_match) && (read_informative_kmer_hits > min_hits_for_informative_read || read_informative_kmer_hits_pe2 > min_hits_for_informative_read)) {
		if (total_read_kmer_hits >= min_hits_for_good_match && total_read_informative_kmer_hits >= min_hits_for_informative_read) {

			/* no longer using the read hits output as was not as useful as the kmer coverage / depth */
			/* printf("%s\t%d\t%d\t%d\t%d\n", PE1_file, read_kmer_hits, read_informative_kmer_hits, read_kmer_hits_pe2, read_informative_kmer_hits_pe2); */

			/* now run through the sequence for loop again and put the read_kmer_hits as the value of *count if count < read_kmer_hits */
			/* only do it when over the threshold old to save time */
			seed_seq = sequence_PE1_copy;
			if (sequence_PE1_length >= seed) { // skip if a read has been trimmed shorter than the seed length
				for (i = 0; i<sequence_PE1_length - seed+1; i++) { // for each possible seed position
					temp_nuc = seed_seq[seed];
					seed_seq[seed] = '\0';
		
					orientStr = orient_string(seed_seq, seedStrRevComp, seed);
					if (!contains_N(orientStr)) {
						count = (unsigned int*)BIO_searchHash(h,orientStr);
						if (count != NULL && count[KMER_TYPE] == INFORMATIVE_KMER) { // only track the stats on the informative kmers
							#ifdef NO_GZIP_OUTPUT
							fprintf(out, "%s\t%d\t%d\t%d\t%d\t%s\n", PE1_file, read_kmer_hits, read_informative_kmer_hits, read_kmer_hits_pe2, read_informative_kmer_hits_pe2, orientStr);
							#else
							gzprintf(gzout, "%s\t%d\t%d\t%d\t%d\t%s\n", PE1_file, read_kmer_hits, read_informative_kmer_hits, read_kmer_hits_pe2, read_informative_kmer_hits_pe2, orientStr);
							#endif
							// TO DO: print individual kmers here in addition to the old way of printing the kmers as a group (so can calculate depth and not just coverage)
							// print it out in same format as below but I think we should print the kmer itself on the end for debugging purposes
							//
							// TO DO: I don't think we need the below lines anymore because they can be calculated from the individual kmer hits
///							if (count[TOTAL_KMER_MATCH_PE1] + count[TOTAL_KMER_MATCH_PE2] < total_read_kmer_hits) // want to keep the max total kmer hit of each kmer, but now doing the read as a unit so one metric cannot creep up while the other stays low; we are optimizing on the total kmer because min informative == 1 performs very well but perhaps optimizing on the sum of PE1 and PE2 will remove junk reads; worse case we go back to the old independent way
//								count[TOTAL_KMER_MATCH_PE1]=read_kmer_hits;
//								count[TOTAL_KMER_MATCH_PE2]=read_kmer_hits_pe2;
//								count[INFORMATIVE_KMER_MATCH_PE1]=read_informative_kmer_hits;
//								count[INFORMATIVE_KMER_MATCH_PE2]=read_informative_kmer_hits_pe2;
	
/* new 	way is trying to look for the best total read across PE1 and PE2 */
/*							if (count[TOTAL_KMER_MATCH_PE1] < read_kmer_hits) // want to keep the max total kmer hit of each kmer
								count[TOTAL_KMER_MATCH_PE1]=read_kmer_hits;
							if (count[INFORMATIVE_KMER_MATCH_PE1] < read_informative_kmer_hits) // want to keep the max informative kmer hit of each kmer
							count[INFORMATIVE_KMER_MATCH_PE1]=read_informative_kmer_hits;
*/
						}
					}
	
					seed_seq[seed] = temp_nuc;
					seed_seq++;
				}
			}


			if (is_PE) {
				seed_seq = seq2->seq.s;
				if (seq2->seq.l >= seed) { // skip if a read has been trimmed shorter than the seed length
					for (i = 0; i<seq2->seq.l - seed+1; i++) { // for each possible seed position
						temp_nuc = seed_seq[seed];
						seed_seq[seed] = '\0';
			
						orientStr = orient_string(seed_seq, seedStrRevComp, seed);
						if (!contains_N(orientStr)) {
							count = (unsigned int*)BIO_searchHash(h,orientStr);
							if (count != NULL && count[KMER_TYPE] == INFORMATIVE_KMER) { // only track the stats on the informative kmers
								#ifdef NO_GZIP_OUTPUT
								fprintf(out, "%s\t%d\t%d\t%d\t%d\t%s\n", PE1_file, read_kmer_hits, read_informative_kmer_hits, read_kmer_hits_pe2, read_informative_kmer_hits_pe2, orientStr);
								#else
								gzprintf(gzout, "%s\t%d\t%d\t%d\t%d\t%s\n", PE1_file, read_kmer_hits, read_informative_kmer_hits, read_kmer_hits_pe2, read_informative_kmer_hits_pe2, orientStr);
								#endif
								// TO DO: same as above for PE1 can remove below and just print out the kmer
								//if (count[TOTAL_KMER_MATCH_PE1] + count[TOTAL_KMER_MATCH_PE2] < total_read_kmer_hits) // want to keep the max total kmer hit of each kmer, but now doing the read as a unit so one metric cannot creep up while the other stays low; we are optimizing on the total kmer because min informative == 1 performs very well but perhaps optimizing on the sum of PE1 and PE2 will remove junk reads; worse case we go back to the old independent way
							//		count[TOTAL_KMER_MATCH_PE1]=read_kmer_hits;
							//		count[TOTAL_KMER_MATCH_PE2]=read_kmer_hits_pe2;
							//		count[INFORMATIVE_KMER_MATCH_PE1]=read_informative_kmer_hits;
							//		count[INFORMATIVE_KMER_MATCH_PE2]=read_informative_kmer_hits_pe2;
							}
						}
		
						seed_seq[seed] = temp_nuc;
						seed_seq++;
					}
				}
			}
		}
	}

	#ifdef NO_GZIP_OUTPUT
	fprintf(out, "#%s\ttotal_kmer_evaluated\t%lld\n", PE1_file, total_kmers_evaluated);
	fprintf(out, "#%s\ttotal_reads_evaluated\t%lld\n", PE1_file, total_reads_evaluated);
	fprintf(out, "#%s\ttotal_genome_kmers\t%lld\n", PE1_file, genome_kmers);
	fprintf(out, "#%s\ttotal_genome_informative_kmers\t%lld\n", PE1_file, genome_informative_kmers);
	#else
	gzprintf(gzout, "#%s\ttotal_kmer_evaluated\t%lld\n", PE1_file, total_kmers_evaluated);
	gzprintf(gzout, "#%s\ttotal_reads_evaluated\t%lld\n", PE1_file, total_reads_evaluated);
	gzprintf(gzout, "#%s\ttotal_genome_kmers\t%lld\n", PE1_file, genome_kmers);
	gzprintf(gzout, "#%s\ttotal_genome_informative_kmers\t%lld\n", PE1_file, genome_informative_kmers);
	#endif

	// TO DO : with the new format can remove this piece of code and just calculate the coverage in python later 
	// still need to do the reset piece when should those be reset? maybe above when we add the new print? need to get that right

	/* done reading the file now print out the kmers that pas the threshold */
/*
	for (i=0; i<hash_size; i++) {
		count = (unsigned int*)BIO_searchHash(h,allKeys[i]);
		// TO DO below is a clear bug that was only look at PE1 AND should have >= and just has > this almost certainly messes up sensitivity
		// but this code will be removed anyways
		if (count[KMER_TYPE] == INFORMATIVE_KMER && count[TOTAL_KMER_MATCH_PE1] > min_hits_for_good_match && count[INFORMATIVE_KMER_MATCH_PE1] > min_hits_for_informative_read) {
			fprintf(out, "%s\t%d\t%d\t%d\t%d\n", PE1_file, count[TOTAL_KMER_MATCH_PE1], count[INFORMATIVE_KMER_MATCH_PE1],count[TOTAL_KMER_MATCH_PE2], count[INFORMATIVE_KMER_MATCH_PE2]);
		}
		// *count = 0; // reset to zero so ready for the next file
		count[TOTAL_KMER_MATCH_PE1] = 0; // reset to zero so ready for the next file
		count[INFORMATIVE_KMER_MATCH_PE1] = 0; // reset to zero so ready for the next file
		count[TOTAL_KMER_MATCH_PE2] = 0; // reset to zero so ready for the next file
		count[INFORMATIVE_KMER_MATCH_PE2] = 0; // reset to zero so ready for the next file
	}
*/

	free(seedStrRevComp);
	free(sequence_PE1_copy);
	gzclose(fp);
	gzclose(fp2);
}


/* name is out of date as this used to hash the informative kmers that passed scrubbing 
   however now we just flag them as type INFORMATIVE_KMER and keep all kmers from the reference genome */
unsigned int hash_scrubbed_kmers(const char *filename, BIO_hash h, const int seed) {
	gzFile fp;
	char *line = NULL;
	size_t max_len = 100;
	unsigned int *count;
	char *pos;
	char *orientStr; // in theory the strings should be oriented but when I get kmers from other software they might not orient the same way
	char *seedStrRevComp = (char*)malloc(sizeof(char) * max_len+1);
	unsigned int num_informative_kmers = 0;
	fp = gzopen(filename, "r");
	line = (char*)malloc(sizeof(char)*max_len);

//	fprintf(stderr, "reading %s to make hash\n", filename);
//	printf("labeling informative kmers\n");
        if (fp == NULL) {
                fprintf(stderr, "could not read file %s in hash_scrubbed_kmers()\n", filename);
                exit(EXIT_FAILURE);
        }

	while (gzgets(fp, line, max_len)) {
		if (line[0] != '#')  { // skip comments
			if ((pos=strchr(line, '\n')) != NULL)
				*pos = '\0';

			// the provided kmers must be the length of the kmers we're searching for
			if (strlen(line) == seed) {
				orientStr = orient_string(line, seedStrRevComp, seed);
//				printf("line\torient\t%s\t%s\t%d\n", line, orientStr, strlen(orientStr));


				// old from when we only hashed the informative kmers
//				count = calloc(vec_size, sizeof(unsigned int));
				
//				count = (unsigned int*)BIO_searchHash(h,line);
				count = (unsigned int*)BIO_searchHash(h,orientStr);
				if (count != NULL) {
//					printf("hit %s\ttype %d\n", line, count[KMER_TYPE]);
					count[KMER_TYPE] = INFORMATIVE_KMER;
					num_informative_kmers++;
//					printf("now type %d\n", count[KMER_TYPE]);
				}
				else {
					printf("error could not find informative kmer %s in the total kmer list\n", line);
				}
			}
			else {
				printf("error string length in the scrubbed kmer file (%s) must be the same size as the kmer length (scrubbed kmer, scrubbed kmer len, seed len): %s, %d, %d\n", filename, line, (int)strlen(line), seed);
			}
		}
	}

//	printf("hash size is %d\n", BIO_getHashSize(h));

	free(line);
	free(seedStrRevComp);
	gzclose(fp);

	return num_informative_kmers;
}

unsigned int get_file_type(const char *filetype) {
//	printf("checking %s\n", filetype);
	if (strcmp(filetype, "SE") == 0 || strcmp(filetype, "se") == 0) {
		return NOT_PAIRED_END;
	}
	else if (strcmp(filetype, "PE") == 0 || strcmp(filetype, "pe") == 0) {
		return IS_PAIRED_END;
	}
	else if (strcmp(filetype, "PEI") == 0 || strcmp(filetype, "pei") == 0 || strcmp(filetype, "IPE") == 0 || strcmp(filetype, "ipe") == 0) {
		return IS_PAIRED_END_INTERLEAVE;
	}
	else {
		return UNKNOWN_FILE_TYPE;
	//	printf("unrecognized file type: %s\n", filetype);
	
	//	exit(0);
	}

	return UNKNOWN_FILE_TYPE;
}


void usage() {
	fprintf(stderr, "Usage paired end with 2 files:\n\tstrain_detect -r <reference_genome.fna> -a <informative_kmer_file.txt> -b <paired-end-file1> -c <paired-end-file1> -t PE -o <kmer outfile>\n");
	fprintf(stderr, "Usage paired end interleaved 1 file:\n\tstrain_detect -r <reference_genome.fna> -a <informative_kmer_file.txt> -b <paired-end-file1>  -t PEI -o <kmer outfile>\n");
	fprintf(stderr, "Usage single end 1 file:\n\tstrain_detect -r <reference_genome.fna> -a <informative_kmer_file.txt> -b <single-end-file1>  -t SE -o <kmer outfile>\n");
	fprintf(stderr, "Usage single end 1 file:\n\tstrain_detect -r <reference_genome.fna> -a <informative_kmer_file.txt> -B <batch-list-of-metagenomes> -o <kmer outfile>\n\n");
	fprintf(stderr, "format for metagenomics batch file is:\n");
	fprintf(stderr, "PE\tfile1_PE1.fasta\tfile1_PE2.fasta\n");
	fprintf(stderr, "SE\tfile1_PE1.fasta\n");
	fprintf(stderr, "PEI\tfile1_PE1.fasta\n");
	fprintf(stderr, "\nlines that begin with # are considered comments and ignored\n");
	fprintf(stderr, "\ninformative kmer file is a list of all of the kmers left in the reference genome post scrubbing\n");
//	fprintf(stderr, "(alternative)         -A <file with multiple filenames> -b <metagenome.fasta|.fastq>\n");
//
	return;
}
