import os
from snakemake.io import expand, directory
import re

# PATH
ALL_GENOMES_DIR = 'strains/'
ALL_METAGENOMES_DIR = 'metagenomes/'
MIN_FRACTION_TO_KEEP = "0.01"
PANGENOME_SCRUB_DB = "genomes_to_scrub.txt"
METAGENOME_SCRUB_DB = "metagenomes_to_scrub.txt"
SCRUB_FILTER = 'python ../scripts/kmer_scrub_filter.py'
MAPPING_TARGETS = "target_metagenomes.txt"
#BACKGROUND_METAGENOMES = 'background_metagenomes.txt'
COVERAGE_SCRIPT = 'python ../scripts/coverage_depth.py'
STRAIN_DETECT = "../src/strain_detect"
MIN_KMER_HITS = '5'
KMER_SCRUB_COUNT = '../src/kmer_scrub_count'


FILES = [i for i in os.listdir(ALL_GENOMES_DIR) if i.endswith('.fna.gz')]
STRAINS = [re.split('.fna', i)[0] for i in sorted(FILES)]

rule all:
    input: expand(ALL_GENOMES_DIR + "{sample}.coverage_depth", sample=STRAINS), MAPPING_TARGETS

rule make_mapping_targets:
    output: MAPPING_TARGETS
    shell: "ls " + ALL_METAGENOMES_DIR + "/* | egrep 'fastq|fna|fastq   > {output}"


rule kmer_scrub_counts:
    input: pan = PANGENOME_SCRUB_DB, meta = METAGENOME_SCRUB_DB, fna = ALL_GENOMES_DIR + "{sample}.fna.gz"
    output: kmer_counts = ALL_GENOMES_DIR + "{sample}.scrub_kmer_counts.gz", progress = ALL_GENOMES_DIR + "{sample}.progress"
    shell: KMER_SCRUB_COUNT +  " -r {input.fna} -A {input.pan} -B {input.meta} -p {output.progress} | gzip --best > {output.kmer_counts}"

rule kmer_scrub:
    input:  ALL_GENOMES_DIR + "{sample}.scrub_kmer_counts.gz"
    output: ALL_GENOMES_DIR + "{sample}.scrubbed_kmers.gz"
    shell: SCRUB_FILTER + " -s {input} --min_fraction " + MIN_FRACTION_TO_KEEP + " | gzip --best > {output}"

rule strain_detect:
    input: kmers = ALL_GENOMES_DIR + "{sample}.scrubbed_kmers.gz", fna = ALL_GENOMES_DIR + "{sample}.fna.gz", targets = MAPPING_TARGETS
    output:  kmer = ALL_GENOMES_DIR + "{sample}.kmer_hits.gz"
    shell: STRAIN_DETECT + " -r {input.fna} -a {input.kmers} -B {input.targets} -o {output.kmer}"

rule coverage_results:
    input: kmer_hits = ALL_GENOMES_DIR + "{sample}.kmer_hits.gz"
    output: ALL_GENOMES_DIR + "{sample}.coverage_depth"
    shell: COVERAGE_SCRIPT + " -k {input.kmer_hits} --min_kmer_hits " + MIN_KMER_HITS + "  > {output}"
#    shell: COVERAGE_SCRIPT + " -k {input.kmer_hits} -b {input.background_meta} --min_kmer_hits " + MIN_KMER_HITS + "  > {output}"

