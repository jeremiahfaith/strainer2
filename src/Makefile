CFLAGS = -O2 -Wall -Wextra -Wshadow -Wcast-qual -g

LDFLAGS = -lz -lm

src = $(wildcard *.c)
obj = $(src:.c=.o)
shared_obj = BIO_sequence.o BIO_hash.o genome_compare.o up2bit.o



#all: genome_compare up2bit_test strain_track
#all: genome_compare strain_track pangenome_compare kmer_matrix kmer_scrub_count strain_detect
all: genome_compare kmer_scrub_count strain_detect

#genome_compare: $(obj)
genome_compare: main.o $(shared_obj)
	$(CC) -o $@ $^ $(LDFLAGS)

kmer_scrub_count: kmer_scrub_count.o $(shared_obj)
	$(CC) -o $@ $^ $(LDFLAGS)

#up2bit_test: up2bit.o
#	$(CC) -o $@ $^ $(LDFLAGS)

strain_detect: strain_detect.o $(shared_obj)
	$(CC) -o $@ $^ $(LDFLAGS)


.PHONY: clean
clean:
	rm $(obj) strain_detect kmer_scrub_count genome_compare
