echo ''
echo 'STEP1: kmer_scrub_count (counting the frequency of target genome kmers in unrelated genomes and metagenomes to learn their frequency'
echo 'example should take 1-3 minutes to run'
../src/kmer_scrub_count -r strains/Bacteroides_ovatus_1001283st1_B8_1001283B150210_160208.fna.gz -A genomes_to_scrub.txt -B metagenomes_to_scrub.txt -p strains/Bacteroides_ovatus_1001283st1_B8_1001283B150210_160208.progress | gzip --best > strains/Bacteroides_ovatus_1001283st1_B8_1001283B150210_160208.scrub_kmer_counts.gz
echo 'STEP1: complete'

echo ''
echo ''
echo 'STEP2: keeping the 1% most rare kmers in the target genome'
echo 'example should take <1 minute to run'
python ../scripts/kmer_scrub_filter.py -s strains/Bacteroides_ovatus_1001283st1_B8_1001283B150210_160208.scrub_kmer_counts.gz -m 0.01 | gzip --best > strains/Bacteroides_ovatus_1001283st1_B8_1001283B150210_160208.scrubbed_kmers.gz
echo 'STEP2: complete'

echo ''
echo ''
echo 'STEP3: identifying informative (rare) kmers in the target metagenomes of interest'
echo 'example should take <<1 minute to run'
../src/strain_detect -r strains/Bacteroides_ovatus_1001283st1_B8_1001283B150210_160208.fna.gz -a strains/Bacteroides_ovatus_1001283st1_B8_1001283B150210_160208.scrubbed_kmers.gz -B target_metagenomes.txt -o strains/Bacteroides_ovatus_1001283st1_B8_1001283B150210_160208.kmer_hits.gz
echo 'STEP3: complete'


echo ''
echo ''
echo 'STEP4: counting the coverage and depth of the informative kmers that were found in each target metagenome'
echo 'example should take <<1 minute to run'
python ../scripts/coverage_depth.py -k strains/Bacteroides_ovatus_1001283st1_B8_1001283B150210_160208.kmer_hits.gz > strains/Bacteroides_ovatus_1001283st1_B8_1001283B150210_160208.coverage_depth
echo 'STEP4: complete'
echo ''
