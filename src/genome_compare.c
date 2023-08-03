#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include "kseq.h"
#include "BIO_sequence.h"
#include "BIO_hash.h"
#include "genome_compare.h"
#include <inttypes.h>
#include <time.h>
#include "up2bit.h"

KSEQ_INIT(gzFile, gzread)

/* 
	TO DO

	1) switch to tommy hash (fixed size) to see if faster

*/

BIO_sequences GEN_read_seq_file(char *file, int *genome_length);
BIO_hash GEN_hash_sequences(BIO_sequences s, const int seed, int genome_length);
inline int rc_strcmp(char *seed_seq, int seed);
//void eliminate_nonunique_keys(BIO_hash seqHash, unsigned int *, unsigned int *);
char* unique_name_suffix(const char *original_name, const char *part1, const char *suffix);
void GEN_hash_background_subtract(const char *subtract_file, BIO_hash seqHash, const int seed, unsigned int *removed_Nmers);
void write_hash_distribution(BIO_hash h, FILE *out);
void write_hash_matrix(BIO_hash h, FILE *out, int num_files, char **file_names);
void sum_and_instances(unsigned int *counts, int count_len, unsigned int *sum, unsigned int *instances);
void GEN_calculate_kmer_count(const char *file, const int seed, BIO_hash h, unsigned int vec_column);

void write_hash_distribution(BIO_hash h, FILE *out) {
	unsigned int i;
	unsigned int *count = NULL;

	for (i=0; i<h->M; i++) {
		count = h->data[i].DATA;
		if (count != NULL) 
			if (*count > 0) 
				fprintf(out, "%d\n", *count);
	}
}

void write_hash_matrix(BIO_hash h, FILE *out, int num_files, char **file_names) {
	unsigned int i,j;
	unsigned int *counts;
	unsigned int min_sum = 4; // needs at least N counts
	unsigned int min_instances = 2; // needs show up in a least N of the possible recipient files
	unsigned int max_instances = 5; // needs show up in a least N of the possible recipient files
	unsigned int sum = 0;
	unsigned int instances = 0;

	fprintf(out, "kmer");
	for (i=0; i<num_files; i++)
		fprintf(out, "\t%s", file_names[i]);
	fprintf(out, "\n");

//	return;
	for (i=0; i<h->M; i++) {
		counts = h->data[i].DATA;
		if (counts != NULL)  {
			sum = 0;
			instances = 0;
			sum_and_instances(counts, num_files, &sum, &instances);
			
			if (sum >= min_sum && instances >= min_instances && instances < max_instances) {
				fprintf(out, "%s", h->data[i].key);
				for (j=0; j<num_files; j++) {
					fprintf(out, "\t%d", counts[j]);
				}
				fprintf(out, "\n");
			}
		}
	}

}

void sum_and_instances(unsigned int *counts, int count_len, unsigned int *sum, unsigned int *instances) {
	unsigned int i;

	for (i=0; i<count_len; i++) {
		if (counts[i])
			*instances += 1;

		*sum += counts[i];
	}
}


void eliminate_nonunique_keys(BIO_hash h, unsigned int *non_unique_count, unsigned int *total_key_count) {
	unsigned int i;
	unsigned int *count = NULL;
	

	for (i=0; i<h->M; i++) {
		count = h->data[i].DATA;
		if (count != NULL) {
			*total_key_count += 1;
			if (*count > 0) {
		//		printf("non-unique %s\t%d\n", h->data[i].key, *count);
				*non_unique_count += 1;

				free(h->data[i].DATA);	
				h->data[i].DATA = NULL;
				free(h->data[i].key);	
	                        h->N--; // decrement the number of keys in hash
//				k->data[i].key = NULL;
			}
		}
	}

}

void GEN_all_kmer_counts_skip_file(const char *B_file, const char *skip_file, const int seed, BIO_hash seqHash, unsigned int vec_column, FILE *progress) {
        FILE *fp;
        char *line = NULL;
        size_t len = 0;
        ssize_t readL;
        fp = fopen(B_file, "r");
        char *pos;
	time_t ltime;

        if (fp == NULL) {
                fprintf(stderr, "could not read file %s in GEN_all_kmer_counts()\n", B_file);
                exit(EXIT_FAILURE);
        }

        while ((readL = getline(&line, &len, fp)) != -1) {
                if ((pos=strchr(line, '\n')) != NULL)
                        *pos = '\0';

		if (progress != NULL) {
			ltime = time(NULL); // current time
			fprintf(progress, "%s\t%s", line, asctime(localtime(&ltime)));
		}

		if (strcmp(skip_file, line) != 0)  // don't hash if it is the genome itself
	                GEN_calculate_kmer_count(line, seed, seqHash, vec_column);
		else 
			fprintf(stderr, "skipping %s (identical match)\n", line);
        }

        fclose(fp);
        free(line);
}


void GEN_all_kmer_counts(const char *B_file, const int seed, BIO_hash seqHash, unsigned int vec_column, FILE *progress) {
        FILE *fp;
        char *line = NULL;
        size_t len = 0;
        ssize_t readL;
        fp = fopen(B_file, "r");
        char *pos;
	time_t ltime;

        if (fp == NULL) {
                fprintf(stderr, "could not read file %s in GEN_all_kmer_counts()\n", B_file);
                exit(EXIT_FAILURE);
        }

        while ((readL = getline(&line, &len, fp)) != -1) {
                if ((pos=strchr(line, '\n')) != NULL)
                        *pos = '\0';

		if (progress != NULL) {
			ltime = time(NULL); // current time
			fprintf(progress, "%s\t%s", line, asctime(localtime(&ltime)));
		}

                GEN_calculate_kmer_count(line, seed, seqHash, vec_column);
        }

        fclose(fp);
        free(line);
}

void GEN_calculate_kmer_count(const char *file, const int seed, BIO_hash h, unsigned int vec_column) {
	gzFile fp;
	kseq_t *seq;
	int l;
	unsigned int i;
	char *seedStrRevComp = (char*)malloc(sizeof(char) * seed+1); 
	char *orientStr;
	unsigned int *count = NULL;
	unsigned int seq_count = 0;
	seedStrRevComp[seed] = '\0';
	char temp_nuc;
	char *seed_seq;
	int enough_reads = 0;

	fp = gzopen(file, "r");
	if (fp == NULL) {
		fprintf(stderr, "could not read file %s in GEN_calculate_kmer_count()\n", file);
		exit(EXIT_FAILURE);
	}
	seq = kseq_init(fp);

//	printf("%s\n", file);

	while (l = kseq_read(seq) >= 0) {
//	while (!enough_reads && (l = kseq_read(seq)) >= 0) {
//		if (max_reads && seq_count++ > max_reads)
//			enough_reads = 1;
		BIO_stringToUpper(seq->seq.s); // keep same case
		seed_seq = seq->seq.s;

		for (i = 0; i<seq->seq.l - seed+1; i++) { // for each possible seed position
			temp_nuc = seed_seq[seed];
			seed_seq[seed] = '\0';	

			orientStr = orient_string(seed_seq, seedStrRevComp, seed);

			if (!contains_N(orientStr)) {
				count = (unsigned int*)BIO_searchHash(h,orientStr);
				if (count != NULL) {
					count[vec_column]+=1; // increment the n-mer
				}

			} 

			seed_seq[seed] = temp_nuc;
			seed_seq++;
		}
	}

	kseq_destroy(seq);
	gzclose(fp);
	free(seedStrRevComp);
}





void GEN_all_coverage(const char *a_file, const char *B_file, const int seed, BIO_hash seqHash, unsigned int max_seeds, double threshold_for_fullmap) {
        FILE *fp;
        char *line = NULL;
        size_t len = 0;
        ssize_t readL;
        fp = fopen(B_file, "r");
        char *pos;

        if (fp == NULL) {
                fprintf(stderr, "could not read file %s in GEN_all_coverage()\n", B_file);
                exit(EXIT_FAILURE);
        }

        while ((readL = getline(&line, &len, fp)) != -1) {
                //printf("%s", line);

                if ((pos=strchr(line, '\n')) != NULL)
                        *pos = '\0';

                GEN_calculate_coverage(a_file, line, seed, seqHash, max_seeds, threshold_for_fullmap);
        }

//      printf("calculating coverage to all in %s\n", B_file);

        fclose(fp);
        free(line);
}


void GEN_calculate_coverage(const char *hash_file, const char *file, const int seed, BIO_hash h, unsigned int max_seeds, double threshold_for_fullmap) {
	gzFile fp;
	kseq_t *seq;
	int l;
	unsigned int i,j;
	int idx;
	int hits = 0;
	int misses = 0;
	char *seedStr        = (char*)malloc(sizeof(char) * seed+1); 
	char *seedStrRevComp = (char*)malloc(sizeof(char) * seed+1); 
	char *orientStr;
	int N_skip;
	int fullmap = 0;
	seedStr[seed] = '\0';
	seedStrRevComp[seed] = '\0';

	fp = gzopen(file, "r");
	if (fp == NULL) {
		fprintf(stderr, "could not read file %s in GEN_calculate_coverage()\n", file);
		exit(EXIT_FAILURE);
	}
	seq = kseq_init(fp);

	while ((l = kseq_read(seq)) >= 0) {
		BIO_stringToUpper(seq->seq.s); // keep same case
//		printf("%s\tlen: %d\n", seq->seq.s, seq->seq.l);
//		continue;
		if (seq->seq.l<seed) // don't compare sequences that are too short
			continue;

		for (i = 0; i<seq->seq.l - seed+1; i++) { // for each possible seed position
			idx=0;
			N_skip = 0;
			for (j=i; j<i+seed; j++) {  // generate the seeds
				seedStrRevComp[idx] = seedStr[idx] = seq->seq.s[j];
				if (seq->seq.s[j] == 'N') // don't hash when we have N
					N_skip = 1;
				idx++;
			}

//			printf("%d len %d\t%d\t%s\n", i, seq->seq.l, seq->seq.l-seed+1, seedStr);

			/* analyze the seed */
			if (!N_skip) {
				BIO_reverseComplement(seedStrRevComp);
				if (strcmp(seedStr, seedStrRevComp) > 0) //  orient the strings
					orientStr = seedStr;
				else
					orientStr = seedStrRevComp;;
	
				if (BIO_searchHash(h,orientStr) != NULL)
					hits++;
				else
					misses++;

			}
			if (max_seeds && hits + misses >= max_seeds) {
//				printf("max hit break %d\n", hits+misses);
//
				if (!fullmap) {
					if (!fullmap && (double)hits/(hits+misses) > threshold_for_fullmap) { // for good scores might want to force the entire score to be more accurate
						fullmap = 1;
				//		printf("doing fullmap\n");
					}
					else { // score is low; move on
						goto end_loop;
						break;
					}
				}
			}
//			printf("%s\t%s\n", seedStr, seedStrRevComp);
		}
	}

	end_loop:;

	printf("%s\t%s\t%d\t%d\t%f\n", hash_file, file, hits, misses, (double)hits/(hits+misses));

//		seq->seq.s; // this is the varible to look at
	kseq_destroy(seq);
	gzclose(fp);
	free(seedStr);
	free(seedStrRevComp);
}

uint64_t GEN_metagenome_coverage_to_ref(const char *file, const int seed, BIO_hash h, uint64_t *num_matches, unsigned int max_reads) {
	gzFile fp;
	kseq_t *seq;
	int l;
	unsigned int i;
//	int i,j;
//	int idx;
//	char *seedStr        = (char*)malloc(sizeof(char) * seed+1); 
	char *seedStrRevComp = (char*)malloc(sizeof(char) * seed+1); 
	char *orientStr;
//	int N_skip;
	unsigned int *count = NULL;
	uint64_t non_N_counts = 0;
//	seedStr[seed] = '\0';
	seedStrRevComp[seed] = '\0';
	unsigned int seq_count = 0;
//	time_t seconds;
//	time_t prev;
	char temp_nuc;
	char *seed_seq;
	int enough_reads = 0;
	//uint64_t up2bit = 0;

	fp = gzopen(file, "r");
	if (fp == NULL) {
		fprintf(stderr, "could not read file %s in GEN_metagenome_coverage_to_ref()\n", file);
		exit(EXIT_FAILURE);
	}
	seq = kseq_init(fp);

//	printf("ref %s to %s\n", hash_file, file);
//	prev = time(NULL);

	while (!enough_reads && (l = kseq_read(seq)) >= 0) {
		if (max_reads && seq_count++ > max_reads)
			enough_reads = 1;
//		if (seq_count++ % 1000 == 0 && seq_count > 3200000) {
/*		if (seq_count++ % 100000 == 0) {
			seconds = time(NULL);
			printf("on line %d time %d\n", seq_count, seconds-prev);
			printf("%s\n", seq->name.s);
			prev = seconds;
		}
*/
		BIO_stringToUpper(seq->seq.s); // keep same case
//		printf("%s\tlen: %d\n", seq->seq.s, seq->seq.l);
		seed_seq = seq->seq.s;

		for (i = 0; i<seq->seq.l - seed+1; i++) { // for each possible seed position
			temp_nuc = seed_seq[seed];
			seed_seq[seed] = '\0';	
//			printf("%s\n", seed_seq);

			orientStr = orient_string(seed_seq, seedStrRevComp, seed);
//			up2bit = encode_DNA_2_bit(orientStr, seed); // encode test here does NOT significantly increase run time (actually ran 3x with and without this and no difference)
//			up2bit -= 20;
//			decode_DNA_2_bit(up2bit, seed, orientStr);


			if (!contains_N(orientStr)) {
			//	printf("%d\n", contains_N(orientStr));
//				count = (unsigned int*)BIO_searchHash64(h,up2bit);
				count = (unsigned int*)BIO_searchHash(h,orientStr);
				if (count != NULL) {
					*count+=1; // increment the n-mer
					*num_matches+=1;
				}

				non_N_counts++; // not really the correct name because I am not looking for N anymore;
			} 
//			else {
//				printf("have N %s\n", orientStr);
//			}

			seed_seq[seed] = temp_nuc;
			seed_seq++;
		}
	}

	kseq_destroy(seq);
	gzclose(fp);
//	free(seedStr);
	free(seedStrRevComp);

	return non_N_counts;
}

int contains_N(char *str) {
	while (*str != '\0') {
		if (*str == 'N')
			return 1;
		str++;
	}

	return 0;
}


BIO_sequences GEN_read_seq_file(char *file, int *genome_length) {
	gzFile fp;
	kseq_t *seq;
	BIO_sequences seqs = BIO_initSequences();
	int l;

	fp = gzopen(file, "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		BIO_seqItem sI = BIO_createSeqItem(seq->name.s,seq->seq.s); // copy the string into proper format
		BIO_stringToUpper(sI->seq); // make everything same case
		*genome_length += sI->seqLength;
		BIO_addSeqItem(seqs, sI);
	}

	// clean up
	kseq_destroy(seq);
	gzclose(fp);
	return seqs;
}

BIO_hash GEN_hash_sequences(BIO_sequences S, const int seed, int genome_len) {
	int initial_size = genome_len * 2; // initialize at 2x the genome size to increase speed!!!!
        BIO_seqItem current;
	BIO_hash h = BIO_initHash(initial_size);
	int i, j;
	char *seedStr        = (char*)malloc(sizeof(char) * seed+1); 
	char *seedStrRevComp = (char*)malloc(sizeof(char) * seed+1); 
	int idx;
	int N_skip = 0;
	seedStr[seed] = '\0';
	seedStrRevComp[seed] = '\0';
	

	for (current=S->seqs; current != NULL; current=current->next) { // for each sequence in a genome (or set of sequences)
//		printf("%s\n", current->seq);
//		printf("%s\t%d\n", current->name, current->seqLength);
		for (i = 0; i<current->seqLength - seed+1; i++) { // for each possible seed position
			idx=0;
		
			N_skip = 0;
			for (j=i; j<i+seed; j++) {  // generate the seeds
				seedStrRevComp[idx] = seedStr[idx] = current->seq[j];
				if (current->seq[j] == 'N') // don't hash when we have N
					N_skip = 1;
				
				idx++;
			}
			
			if (!N_skip) {
				BIO_reverseComplement(seedStrRevComp);
				if (strcmp(seedStr, seedStrRevComp) > 0) { //  orient the strings
					if (BIO_searchHash(h,seedStr) == NULL) // only add if not there already
						BIO_addHashKeyOnly(h, seedStr);
				}
				else {
					if (BIO_searchHash(h,seedStrRevComp) == NULL)
						BIO_addHashKeyOnly(h, seedStrRevComp);
				}
				
			}
		}
	}
	
	free(seedStr);
	free(seedStrRevComp);
	return h;
}

//void GEN_print_coverage_to_ref(BIO_sequences S, const int seed, BIO_hash h, uint64_t seed_count) {
void GEN_print_coverage_to_ref(const char *file, const int seed, BIO_hash h, FILE *out, unsigned int *used_seeds, unsigned int *possible_seeds, uint64_t *total_counts, const int print_track) {
	gzFile fp;
	kseq_t *seq;
        //BIO_seqItem current;
	//int i, j;
	unsigned int i;
	char *seedStr        = (char*)malloc(sizeof(char) * seed+1); 
	char *seedStrRevComp = (char*)malloc(sizeof(char) * seed+1); 
	//int idx;
//	int N_skip = 0;
	int l;
	seedStr[seed] = '\0';
	seedStrRevComp[seed] = '\0';
	unsigned int *count;
	char temp_nuc;
	char *seed_seq;
	char *orientStr;
	//uint64_t total_counts = 0;
	//unsigned int used_seeds = 0;
	//unsigned int possible_seeds = 0;
	*total_counts = 0;
	*used_seeds = 0;
	*possible_seeds = 0;

	

	fp = gzopen(file, "r");
	if (fp == NULL) {
		fprintf(stderr, "could not read file: %s in GEN_print_coverage_to_ref()\n", file);
		exit(EXIT_FAILURE);
	}
	seq = kseq_init(fp);

//	for (current=S->seqs; current != NULL; current=current->next) { // for each sequence in a genome (or set of sequences)
        while ((l = kseq_read(seq)) >= 0) {
		BIO_stringToUpper(seq->seq.s); // keep same case
		seed_seq = seq->seq.s;

		for (i = 0; i<seq->seq.l - seed+1; i++) { // for each possible seed position
			temp_nuc = seed_seq[seed];
			seed_seq[seed] = '\0';  
			orientStr = orient_string(seed_seq, seedStrRevComp, seed);
			if (!contains_N(orientStr)) {
				count = (unsigned int*)BIO_searchHash(h,orientStr);
				if (count != NULL) {
					if (print_track)
						fprintf(out, "%s\t%d\n", orientStr, *count);

					*total_counts += *count;
					*possible_seeds += 1; // all seeds without an "N"

					if (*count > 0)
						*used_seeds += 1;
				}
				else {
					if (print_track)
						fprintf(out, "%s\t-1\n", orientStr); // this should not happen except on shared eliminated seeds; print -1? 
				}
			}
			else {
				if (print_track)
					fprintf(out, "%s\t-2\n", orientStr); // print -2 for seed with an N
			}
			seed_seq[seed] = temp_nuc;
			seed_seq++;

		}
	}

	if (print_track)
		fprintf(out, "#total_counts\t%" PRIu64 "\n", *total_counts);


	free(seedStr);
	free(seedStrRevComp);
}
void GEN_hash_all_sequences_kmer_mat(const char *A_file,  const int seed) {
        FILE *fp;
        char *line = NULL;
        size_t len = 0;
        ssize_t read_s;
	//FILE *fout;
	char *pos;
	//char *outfile;
	char **infiles = NULL;
	unsigned int num_infiles=0;
	int i;
//	char *suffix = "pangenome";
//	char *suffixB = "pangenome_dist";
	BIO_hash seqHash = BIO_initHash(DEFAULT_GENOME_HASH_SIZE);

        fp = fopen(A_file, "r");
        if (fp == NULL) {
                fprintf(stderr, "could not read file %s in GEN_hash_all_sequences_kmer_mat()\n", A_file);
                exit(EXIT_FAILURE);
        }

	/* count the number of sequence input files */
        while ((read_s = getline(&line, &len, fp)) != -1) 
		num_infiles++;
	rewind(fp);
	infiles = malloc(num_infiles * sizeof(char*));
	
	i=0;
        while ((read_s = getline(&line, &len, fp)) != -1) {
		if ((pos=strchr(line, '\n')) != NULL)
			*pos = '\0';

		fprintf(stderr, "reading file %s\t%d of %d\n", line, i+1, num_infiles);
//		printf("reading file %s\t%d of %d\n", line, i+1, num_infiles);
		GEN_hash_sequences_set_count_vec(line, seed, seqHash, 1, 1, i, num_infiles);
		infiles[i++] = strdup(line); // save the infiles; these will be the headers
        }

	write_hash_matrix(seqHash, stdout, num_infiles, infiles);



	for (i=0; i<num_infiles; i++)
		free(infiles[i]);
//		printf("%s\n", infiles[i]);


	free(infiles);
}

//void GEN_hash_all_sequences_pangenome(const char *A_file, const char *b_file, const int seed) {
void GEN_hash_all_sequences_pangenome(const char *A_file,  const int seed, const char *ref_file, int write_dist) {
        FILE *fp;
        char *line = NULL;
        size_t len = 0;
        ssize_t read_s;
	FILE *fout;
	char *pos;
	char *outfile;
	char *suffix = "pangenome";
	char *suffixB = "pangenome_dist";
	BIO_hash seqHash = BIO_initHash(DEFAULT_GENOME_HASH_SIZE);
	unsigned int used_seeds; 
	unsigned int possible_seeds; 
	uint64_t total_counts;
	int hashed_files = 0;

        fp = fopen(A_file, "r");
        if (fp == NULL) {
                fprintf(stderr, "could not read file %s GEN_hash_all_sequences_pangenome()\n", A_file);
                exit(EXIT_FAILURE);
        }

	/* hash all of the sequences */
        while ((read_s = getline(&line, &len, fp)) != -1) {

                if ((pos=strchr(line, '\n')) != NULL)
                        *pos = '\0';

		fprintf(stderr, "hashing %s\n", line);
		GEN_hash_sequences_set_count(line, seed, seqHash, 1, 1);
		hashed_files++;
//		GEN_hash_sequences_set_count(line, seed, seqHash, 1, 0);
        }

	rewind(fp);
//	printf("hereA\n");
	
	if (ref_file) { // just make track to the one requested
		outfile        = unique_name_suffix(ref_file, "", suffix); 
		printf("file %s to %s\n", ref_file, outfile);
        	fout = fopen(outfile, "w");
        	if (fout == NULL) {
        	        fprintf(stderr, "could not open ref_file: %s in GEN_hash_all_sequences_pangenome()\n", outfile);
        	        exit(EXIT_FAILURE);
		}
		fprintf(fout, "#%s\n", ref_file);
		fprintf(fout, "#output to %s\n", outfile);
		fprintf(fout, "#pangenome_size\t%d\n", hashed_files);
		GEN_print_coverage_to_ref(ref_file, seed, seqHash, fout, &used_seeds, &possible_seeds, &total_counts, 1);
		fclose(fout);
		free(outfile);
	}
	else {
	        while ((read_s = getline(&line, &len, fp)) != -1) { // default to making all tracks
			if ((pos=strchr(line, '\n')) != NULL)
				*pos = '\0';


			outfile        = unique_name_suffix(line, "", suffix); 
			printf("file %s to %s\n", line, outfile);
	        	fout = fopen(outfile, "w");
	        	if (fout == NULL) {
	        	        fprintf(stderr, "could not read file outfile %s in GEN_hash_all_sequences_pangenome()\n", outfile);
	        	        exit(EXIT_FAILURE);
			}
			fprintf(fout, "#%s\n", line);
			fprintf(fout, "#output to %s\n", outfile);
			fprintf(fout, "#pangenome_size\t%d\n", hashed_files);
	//		GEN_print_coverage_to_ref(line, b_file, seed, seqHash, metagenome_seed_count, fout);
			GEN_print_coverage_to_ref(line, seed, seqHash, fout, &used_seeds, &possible_seeds, &total_counts, 1);
			fclose(fout);
			free(outfile);
		}
	}

//	printf("hereB\n");

	if (write_dist) {
		outfile        = unique_name_suffix(A_file, "", suffixB); 
		printf("writing dist to %s\n", outfile);
        	fout = fopen(outfile, "w");
        	if (fout == NULL) {
        	        fprintf(stderr, "could not read file outfile (second) %s in GEN_hash_all_sequences_pangenome()\n", outfile);
        	        exit(EXIT_FAILURE);
		}
		write_hash_distribution(seqHash, fout);
		free(outfile);
		fclose(fout);
	}

        fclose(fp);
	BIO_destroyHash(seqHash); 
        free(line);
}


void GEN_hash_all_sequences_set_count_metagenomics(const char *A_file, const char *b_file, const int seed, const int print_track, unsigned int max_reads) {
 //       BIO_sequences seqs;
        FILE *fp;
        char *line = NULL;
        size_t len = 0;
        ssize_t readL;
	char *pos;
	//int genome_length = 0;
	uint64_t metagenome_seed_count;
	FILE *fout;
	char *outfile;
	BIO_hash seqHash = BIO_initHash(DEFAULT_GENOME_HASH_SIZE);
	unsigned int non_unique_count = 0;
	unsigned int total_key_count = 0;
	unsigned int used_seeds; 
	unsigned int possible_seeds; 
	uint64_t total_counts;
	uint64_t num_matches = 0;
	int num_strains = 0;
	GEN_strainResult *SR = NULL;
	int i=0;
	double scale_sum = 0.0;

        fp = fopen(A_file, "r");
        if (fp == NULL) {
                fprintf(stderr, "could not read file A_file %s in GEN_hash_all_sequences_set_count_metagenomics()\n", A_file);
                exit(EXIT_FAILURE);
        }

	/* hash all of the sequences */
        while ((readL = getline(&line, &len, fp)) != -1) {
              //  printf("%s", line);

                if ((pos=strchr(line, '\n')) != NULL)
                        *pos = '\0';

//		printf("hashing %s\n", line);
		
//		seqs = GEN_read_seq_file(line, &genome_length);
		GEN_hash_sequences_set_count(line, seed, seqHash, 0, 1);
		num_strains++;
		//seqHash = GEN_hash_sequences_set_count(seqs, seed, genome_length);
        }

	SR = malloc(num_strains * sizeof(GEN_strainResult));
	for (i=0; i<num_strains; i++)  {
		SR[i].used_seeds = SR[i].possible_seeds = SR[i].total_counts = 0;
		SR[i].strain_name = NULL;
	}

	eliminate_nonunique_keys(seqHash, &non_unique_count, &total_key_count);

	fprintf(stderr, "eliminate nonunique %d of %d (%f)\n", non_unique_count, total_key_count, (double)non_unique_count/total_key_count);


	/* count metagenome against hash */
	metagenome_seed_count = GEN_metagenome_coverage_to_ref(b_file, seed, seqHash, &num_matches, max_reads);

	/* print counts to each reference */
	rewind(fp);
	i=0;
        while ((readL = getline(&line, &len, fp)) != -1) {
                if ((pos=strchr(line, '\n')) != NULL)
                        *pos = '\0';

//		printf("coverage %s\n", line);
//		printf("output to %s\n", outfile);

		if (print_track) {
			outfile        = unique_name_suffix(line, b_file, "strain_track"); 
	        	fout = fopen(outfile, "w");
	        	if (fout == NULL) {
	        	        fprintf(stderr, "could not read file outfile %s in GEN_hash_all_sequences_set_count_metagenomics()\n", outfile);
	        	        exit(EXIT_FAILURE);
			}

			fprintf(stderr, "output to %s\n", outfile);
			GEN_print_coverage_to_ref(line, seed, seqHash, fout, &used_seeds, &possible_seeds, &total_counts, print_track);
			fprintf(fout, "%s\n", line);
//			fprintf(fout, "output to %s\n", outfile);
			fprintf(fout, "%" PRIu64 "\n", metagenome_seed_count);
			fclose(fout);
			free(outfile);
		}
		else {
			GEN_print_coverage_to_ref(line, seed, seqHash, stderr, &used_seeds, &possible_seeds, &total_counts, print_track);
//			fprintf(stderr, "no print track\n");
		}
		
//		                printf("%s\t%s\t%d\t%d\t%" PRIu64 "\t%" PRIu64 "\t%f\t%f\t%f\n", a_file, b_file, used_seeds, possible_seeds, total_counts, metagenome_seed_count, (double)used_seeds/possible_seeds, (double)total_counts/                  metagenome_seed_count, (double)total_counts/num_matches);

//		printf("%s\t%s\t%d\t%d\t%" PRIu64 "\t%" PRIu64 "\t%f\t%f\t%f\n", line, b_file, used_seeds, possible_seeds, total_counts, metagenome_seed_count, (double)used_seeds/possible_seeds, (double)total_counts/metagenome_seed_count, (double)total_counts/num_matches);
		SR[i].used_seeds = used_seeds;
		SR[i].possible_seeds = possible_seeds;
		SR[i].total_counts = total_counts;
		SR[i].strain_name = strdup(line);
		scale_sum += (double)SR[i].total_counts / SR[i].possible_seeds; // scale by the number of unique seeds in the genome (normalizes by informative genome content so larger genomes do not get higher fractional abundances 
//		printf("%f to scale sum %f\n", (double)SR[i].total_counts/SR[i].possible_seeds, scale_sum);
		i++;
	}

	/* print out the results */
	printf("#query\ttarget\tused_seeds\tpossible_seeds\tseed_counts\tmetagenomic_counts\tfrac_used_seeds\tfrac_counts\tfrac_matches\tscaled_matches\n");
	for (i=0; i<num_strains; i++) {
		printf("%s\t%s\t%d\t%d\t%" PRIu64 "\t%" PRIu64 "\t%f\t%f\t%f\t%f\n", SR[i].strain_name, b_file, SR[i].used_seeds, SR[i].possible_seeds, SR[i].total_counts, metagenome_seed_count, (double)SR[i].used_seeds/possible_seeds, (double)SR[i].total_counts/metagenome_seed_count, (double)SR[i].total_counts/num_matches, ((double)SR[i].total_counts/SR[i].possible_seeds)/scale_sum);
	}




        fclose(fp);
	BIO_destroyHash(seqHash); 
        free(line);

	for (i=0; i<num_strains; i++)
		free(SR[i].strain_name);
	free(SR);
}

// the string from this must be freed! 
char* unique_name_suffix(const char *original_name, const char *part1, const char *suffix) {
	char *new_name = (char*)malloc(sizeof(char) * (strlen(original_name) + strlen(part1) + strlen(suffix) + 3)); 
	strcpy(new_name, original_name);
	strcat(new_name, "_");
	strcat(new_name, part1);
	strcat(new_name, ".");
	strcat(new_name, suffix);

	return new_name;
}

void GEN_hash_background_subtract_all(const char *S_file, BIO_hash seqHash, const int seed, unsigned int *removed_Nmers) {
        FILE *fp;
        char *line = NULL;
        size_t len = 0;
        ssize_t readL;
	char *pos;
	//int genome_length = 0;
	uint64_t metagenome_seed_count;


        fp = fopen(S_file, "r");
        if (fp == NULL) {
                fprintf(stderr, "could not read file %s in GEN_hash_background_subtract_all()\n", S_file);
                exit(EXIT_FAILURE);
        }

	/* hash all of the sequences */
        while ((readL = getline(&line, &len, fp)) != -1) {
              //  printf("%s", line);

                if ((pos=strchr(line, '\n')) != NULL)
                        *pos = '\0';

		printf("subtracting %s\n", line);
		
        }

	fclose(fp);
}

/* this current version is an absolute cut; if a seed is found once it is eliminated
   a better cut might be to screen N seeds in M files (e.g., N = 1000000 and M = 100) and remove seeds
   that show up in more than X donors or more than Y times
*/
void GEN_hash_background_subtract(const char *subtract_file, BIO_hash h, const int seed, unsigned int *removed_Nmers) {
	unsigned int i;
	char *seedStr        = (char*)malloc(sizeof(char) * seed+1); 
	char *seedStrRevComp = (char*)malloc(sizeof(char) * seed+1); 
	char *orientStr;
        char *seed_seq;
        char temp_nuc;
	int idx;
	int l;
	//int N_skip = 0;
	unsigned int *count = NULL;
        gzFile fp;
        kseq_t *seq;

	seedStr[seed] = '\0';
	seedStrRevComp[seed] = '\0';

	fp = gzopen(subtract_file, "r");
	if (fp == NULL) {
		fprintf(stderr, "could not read background subtract file %s in GEN_hash_background_subtract()\n", subtract_file);
		exit(EXIT_FAILURE);
	}
	seq = kseq_init(fp);
	
        while ((l = kseq_read(seq)) >= 0) {
                BIO_stringToUpper(seq->seq.s); // keep same case
                seed_seq = seq->seq.s;

		for (i = 0; i<seq->seq.l - seed+1; i++) { // for each possible seed position
			idx=0;
			temp_nuc = seed_seq[seed];
			seed_seq[seed] = '\0';  
			orientStr = orient_string(seed_seq, seedStrRevComp, seed);

			if (!contains_N(orientStr)) {
				count = (unsigned int*)BIO_searchHash(h,orientStr);
				if (count == NULL) {// only add if not there already
//					BIO_addHashData(h, orientStr, count);
					BIO_removeHashKey(h, orientStr);
					free(count); // free up the memory of the data in the hash
					*removed_Nmers += 1;
				}
			}
		
			seed_seq[seed] = temp_nuc;
			seed_seq++;
		}
	}

	free(seedStr);
	free(seedStrRevComp);
        kseq_destroy(seq);
        gzclose(fp);
}

void GEN_hash_sequences_set_count_vec(const char *file, const int seed, BIO_hash h, const int default_count, const int increment, const int vec_idx, const int vec_size) {
	unsigned int i;
	char *seedStr        = (char*)malloc(sizeof(char) * seed+1); 
	char *seedStrRevComp = (char*)malloc(sizeof(char) * seed+1); 
	char *orientStr;
        char *seed_seq;
        char temp_nuc;
	int idx;
	int l;
	//int N_skip = 0;
	unsigned int *counts = NULL;
        gzFile fp;
        kseq_t *seq;

	seedStr[seed] = '\0';
	seedStrRevComp[seed] = '\0';

	fp = gzopen(file, "r");
	if (fp == NULL) {
		fprintf(stderr, "could not read file %s GEN_hash_sequences_set_count_vec()\n", file);
		exit(EXIT_FAILURE);
	}
	seq = kseq_init(fp);

//	printf("vec idx %d\n", vec_idx);
	


        while ((l = kseq_read(seq)) >= 0) {
                BIO_stringToUpper(seq->seq.s); // keep same case
                seed_seq = seq->seq.s;

	//	for (i = 0; i<current->seqLength - seed+1; i++) { // for each possible seed position
		for (i = 0; i<seq->seq.l - seed+1; i++) { // for each possible seed position

			idx=0;
			temp_nuc = seed_seq[seed];
			seed_seq[seed] = '\0';  
			orientStr = orient_string(seed_seq, seedStrRevComp, seed);

			if (!contains_N(orientStr)) {
				counts = (unsigned int*)BIO_searchHash(h,orientStr);
				if (counts == NULL) {// only add if not there already
//					printf("adding %s %d size %d\n", orientStr, vec_idx, vec_size);
					counts        = calloc(vec_size, sizeof(unsigned int));
					counts[vec_idx] = default_count;
					BIO_addHashData(h, orientStr, counts);
				}
				else {
					counts[vec_idx] += increment;
//					printf("increment %s %d size %d count %d\n", orientStr, vec_idx, vec_size, counts[vec_idx]);
				}
			}
		
			seed_seq[seed] = temp_nuc;
			seed_seq++;
		}
	}

	free(seedStr);
	free(seedStrRevComp);
        kseq_destroy(seq);
        gzclose(fp);
}



//BIO_hash GEN_hash_sequences_set_count(BIO_sequences S, const int seed, int genome_len) {
void GEN_hash_sequences_set_count(const char *file, const int seed, BIO_hash h, const int default_count, const int increment) {
//	int initial_size = genome_len * 2; // initialize at 2x the genome size to increase speed!!!!
        //BIO_seqItem current;
//	BIO_hash h = BIO_initHash(initial_size);
	//int i, j;
	unsigned int i;
	char *seedStr        = (char*)malloc(sizeof(char) * seed+1); 
	char *seedStrRevComp = (char*)malloc(sizeof(char) * seed+1); 
	char *orientStr;
        char *seed_seq;
        char temp_nuc;
	int idx;
	int l;
	//int N_skip = 0;
	unsigned int *count = NULL;
        gzFile fp;
        kseq_t *seq;

	seedStr[seed] = '\0';
	seedStrRevComp[seed] = '\0';

	fp = gzopen(file, "r");
	if (fp == NULL) {
		fprintf(stderr, "could not read file %s in GEN_hash_sequences_set_count()\n", file);
		exit(EXIT_FAILURE);
	}
	seq = kseq_init(fp);
	

//	for (current=S->seqs; current != NULL; current=current->next) { // for each sequence in a genome (or set of sequences)
        while ((l = kseq_read(seq)) >= 0) {
                BIO_stringToUpper(seq->seq.s); // keep same case
                seed_seq = seq->seq.s;

	//	for (i = 0; i<current->seqLength - seed+1; i++) { // for each possible seed position
		for (i = 0; i<seq->seq.l - seed+1; i++) { // for each possible seed position

			idx=0;
			temp_nuc = seed_seq[seed];
			seed_seq[seed] = '\0';  
			orientStr = orient_string(seed_seq, seedStrRevComp, seed);

			if (!contains_N(orientStr)) {
				count = (unsigned int*)BIO_searchHash(h,orientStr);
				if (count == NULL) {// only add if not there already
					count        = (unsigned int*)malloc(sizeof(unsigned int));
					*count = default_count;
					BIO_addHashData(h, orientStr, count);
				}
				else if (increment) {
					*count += 1; /* use this to identify duplicate regions and regions occuring in more than one genome */
				}
			}
		
			seed_seq[seed] = temp_nuc;
			seed_seq++;
		}
	}

	free(seedStr);
	free(seedStrRevComp);
        kseq_destroy(seq);
        gzclose(fp);
}

char* orient_string(char *seed_seq, char *seedStrRevComp, int seed) {
	//char *orient;
	int i;
	const char *seedPtr = seed_seq;
	
//	printf("\n\n%s\n", seed_seq);
//	printf("%d\n", rc_strcmp(seed_seq, seed));

	if (rc_strcmp(seed_seq, seed) >= 0) { /* forward sequence wins */
		return seed_seq; // no need to get reverse complement
	}
	else { /* build the reverse complement */
		seedStrRevComp[seed] = '\0';
		for (i=seed-1; i>=0; i--) {
			seedStrRevComp[i] = COMPLEMENT[*seedPtr];
			seedPtr++;
		}
//		printf("%s\trc\n", seedStrRevComp);
		return seedStrRevComp;
	}
}

int rc_strcmp(char *seed_seq, int seed) {
	int i;
	int revIdx;
	char rc;
	char *seedPtr = seed_seq;

	for (i=0; i<seed; i++) {
		revIdx = seed-(i+1);
	
		rc = COMPLEMENT[seed_seq[revIdx]];
//		printf("%c%c\t%d%d\t", *seedPtr, rc, *seedPtr, rc);
		if (*seedPtr > rc)
			return 1;
		else if (rc > *seedPtr)
			return -1;

		seedPtr++;
	}
	return 0; // palindrome
}
