#ifndef GENOME_COMPARE_H
#define GENOME_COMPARE_H 1

/** 
 \file  genome_compare.h
 \brief genome hashing utilities
 
 genome_compare.h contains basic functions for hashing DNA
 

 */

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>


#define DEFAULT_GENOME_HASH_SIZE 8000000

typedef struct strainRes {
	unsigned int used_seeds;
	unsigned int possible_seeds;
	uint64_t total_counts;
	char *strain_name;
} GEN_strainResult;


inline int contains_N(char *str);
inline char* orient_string(char *seed_seq, char *seedStrRevComp, int seed);
BIO_sequences GEN_read_seq_file(char *file, int *genome_length);
BIO_hash GEN_hash_sequences(BIO_sequences s, const int seed, int genome_length);
void GEN_calculate_coverage(const char *a_file, const char *file, const int seed, BIO_hash h, unsigned int max_seeds, double threshold_for_fullmap);
void GEN_all_coverage(const char *a_file, const char *B_file, const int seed, BIO_hash seqHash, unsigned int max_seeds, double threshold_for_fullmap);
//uint64_t GEN_metagenome_coverage_to_ref(const char *file, const int seed, BIO_hash h);
uint64_t GEN_metagenome_coverage_to_ref(const char *file, const int seed, BIO_hash h, uint64_t *num_hits, unsigned int max_reads);
//BIO_hash GEN_hash_sequences_set_count(BIO_sequences S, const int seed, int genome_len);
void GEN_hash_sequences_set_count(const char *file, const int seed, BIO_hash seqHash, const int default_val, const int increment);
void GEN_hash_sequences_set_count_vec(const char *file, const int seed, BIO_hash seqHash, const int default_val, const int increment, const int idx, const int size);
void GEN_hash_all_sequences_set_count_metagenomics(const char *A_file, const char *b_file, const int seed, const int print_track, unsigned int max_reads);
//void GEN_hash_all_sequences_pangenome(const char *A_file, const char *b_file, const int seed);
void GEN_hash_all_sequences_pangenome(const char *A_file, const int seed, const char *ref_to, int write_dist);
void GEN_hash_all_sequences_kmer_mat(const char *A_file, const int seed);
void GEN_hash_all_metagenome_track(const char *A_file, const char *b_file, const int seed);
void GEN_hash_background_subtract_all(const char *B_file, BIO_hash seqHash, const int seed, unsigned int *removed_Nmers);

/* use this function to take an original hash and see if other fast files have similar kmers (and count them) */
void GEN_all_kmer_counts(const char *B_file, const int seed, BIO_hash seqHash, unsigned int vec_column, FILE *progress);

/* similar to above function but allows you to specify a file to skip */
void GEN_all_kmer_counts_skip_file(const char *B_file, const char *skip_file, const int seed, BIO_hash seqHash, unsigned int vec_column, FILE *progress);




void eliminate_nonunique_keys(BIO_hash seqHash, unsigned int *, unsigned int *);




//void GEN_print_coverage_to_ref(BIO_sequences s, const int seed, BIO_hash seqHash, uint64_t seed_count);
//void GEN_print_coverage_to_ref(const char *s, const char *t, const int seed, BIO_hash seqHash, uint64_t seed_count, FILE *out);
//void GEN_print_coverage_to_ref(const char *file, const char *target, const int seed, BIO_hash h, FILE *out, unsigned int *used_seeds, unsigned int *possible_seeds, uint64_t *total_counts); 
void GEN_print_coverage_to_ref(const char *file, const int seed, BIO_hash h, FILE *out, unsigned int *used_seeds, unsigned int *possible_seeds, uint64_t *total_counts, const int print_track); 






#ifdef __cplusplus
}
#endif

#endif // GENOME_COMPARE_H
