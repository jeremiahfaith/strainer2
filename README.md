# strainer2
kmer-based software for detecting bacterial strain genomes inside metagenomes

The Strainer2 software is an updated version of the [Strainer workflow](https://bitbucket.org/faithj02/strainer-metagenomics/src/master/README.md)
that identifies rare kmers in a bacterial strain genome and searches for those rare kmers in one or more metagenomes to determine if a strain
is present or absent in a metagenome.

The program can be used to track bacterial strains from defined live biotherapeutic products, sequenced culture isolates from fecal transplant donors and recipients, or (less tested) MAGs generated from metagenomes of microbial therapeutic recipients or other individuals that might have strain overlap.

If you use Strainer2, please cite:
[**Precise quantification of bacterial strains after fecal microbiota transplantation delineates long-term engraftment and explains outcomes.**](https://doi.org/10.1038/s41564-021-00966-0). Varun Aggarwala, Ilaria Mogno, Zhihua Li, Chao Yang, Graham Britton, Alice Chen-Liaw, Josephine Mitcham, Gerold Bongers, Dirk Gevers, Jose Clemente, Jean-Frederic Colombel, Ari Grinspan, and Jeremiah Faith. Nature Microbiology 2021.

## What's new in version 2
* computationally intensive programs are written in C
* all software features can be controlled through batch files or command line parameters to facilitate automation (e.g., snakemake)


## Installation
in the source directory is a simple Makefile. 

### External libraries and C headers used
* zlib (must be installed on your system, but very likely it is)
* kseq.h (excellent fasta/fastq parser distributed with this software) [original site]()

### How to install
From inside the src directory type:
```
make
```

and it should generate the executables: kmer_scrub_count, strain_detect, and (not necessary for Strainer2 but a helpful tool) genome_compare.

There are a two python files in the scripts directory that are also needed.


## Overview of the algorithm
1. for all kmers in a strain, count how often each kmer occurs in a large set of genomes and metagenomes where the strain is not expected (i.e., unrelated individuals)
2. use the kmer counts from above to keep only the rare kmers and remove the other kmers. the more you scrub the lower the sensitivity and higher the precision. The default (empirically determined) is to keep 1% of the kmers
3. take the scrubbed kmer set and track the kmers in the metagenomes
4. determine if a strain is present or absent based on the frequency of the rare kmers in the metagenome

Note that the software examples and applications are currently set up for differentiating human gut strains in the human gut metagenome. The main ideas should apply to other body sites, but the learning of informative kmers would require metagenomes and genomes from the target sites.

## Description of programs and order of execution

### timing
Each strain is independently tracked, so the sequential steps for each strain genome can be run in parallel to each other. The number of parallel jobs will be the number of strains to be tracked. However, there are many other opportunities for splitting up further if desired. The current version finishes in 1-4 days depending on the number of metagenomes. Most of the time is spent during the initial kmer analysis or in the strain_detect step if there are many metagenomes. If the work flow was broken down not by strain but by strain and by metagenome (and the pieces combined back together at the end), the run time should run in a few hours instead of a few days. Lastly, the examples are watered down to finish in a few minutes. However, they are not sufficient for actual strain tracking (you would have many false positives if you do not study enough metagenomes/genomes to properly rank rare from common kmers).


### programs
Each program if run with no parameters will provide a Usage statement of all the parameters available.


* `kmer_scrub_count`
	* input: takes a strain genome (the one to be tracked) in fasta format as a parameter -r, a file with a list of metagenomes -B, and a file with a list of genomes -A (all in fasta format with or without gzip compression) and a name of a progress file -p. The program counts the frequency of the strain genome kmers inside the set of genomes and metagenomes. It updates the progress file when it completes a genome (to track progress) and prints the final kmer counts separately for the metagenome and genome. Optionally the program will accept an additional -C file that contains a set of genomes that are also strains from the same drug or FMT donor. These co-occuring strains can be quantified separately, as it might be useful to not have overlapping informative kmers used to track multiple strains in the same drug or FMT.
	* output: the complete list of the strains' kmers and the frequency of each kmer in the genome list, metagenome list, and (optionally) drug/FMT strain list.
* `kmer_scrub_filter.py`
	* input: the output file of `kmer_scrub_count` 
	* output: a file of informative (i.e., scrubbed) kmers using the threshold of --min_fraction the parameter --independent will scrub the pangenome and metagenome independently for a more stringent criterion; however this can lead to very few kmers if the intersection between these two is small. The joint scrubbing is the default.
	* joint scrubbing algorithm: to make sure differences in the number of genomes/metagenomes scrubbed aren't a main influence on the kmers that are scrubbed the kmers for each are converted into the frequencies (i.e., metagenome_kmer_count / sum_metagenome_kmer_counts; genome_kmer_count / sum_genome_kmer_count); then we go through all kmers in the genome and rank them by their MAX frequency between metagenome and genome (note that metagenome frequencies are more strongly distributed so they get they dominate the first kmers but then do not contribute much after that because the genome kmers are more evenly distributed; some interesting biology there...); the kmers are then removed from most frequent to least frequent until there are --min_fraction remaining
* `strain_detect`
        * input: takes a reference genome file -r, the informative kmer file from ``kmer_scrub_filter`, and -B file with multiple metagenomics files. can handle paired-end, single-end, or paired end interleaved. you need to specify the file type in the batch file or you can run individually through the commandline.
        * output: a table where each line is an individual informative kmer found in a specific metagenome
* `coverage_depth.py`
	* input: the informative kmer matches from strain_detect
        * output: the proportion of a strains informative kmers covered and the average depth at which the kmers are covered which can be used to decide if a strain is present or absent in a sample




## Documentation to do list
* add a simple example with and without Snakemake that runs in <5 minutes

