echo 'TEST1: strain_detect paired-end'
#time ../src/strain_detect -r Bacteroides_ovatus_1001283st1_B8_1001283B150210_160208.fna.gz -a Bacteroides_ovatus_1001283st1_B8_1001283B150210_160208.scrubbed_kmers.gz -b 1001283B150225_150804_H5_s07_small_PE1.fasta.gz -c 1001283B150225_150804_H5_s07_small_PE2.fasta.gz -t PE -o outfile.strain_detect.test.gz
echo 'running diff to compare strain_detect paired-end output with expected result'
diff outfile.strain_detect.gz outfile.strain_detect.test.gz

echo ''
echo ''
echo 'TEST2 strain_detect single-end'
#time ../src/strain_detect -r Bacteroides_ovatus_1001283st1_B8_1001283B150210_160208.fna.gz -a Bacteroides_ovatus_1001283st1_B8_1001283B150210_160208.scrubbed_kmers.gz -b 1001283B150225_150804_H5_s07_small_PE1.fasta.gz -t SE -o outfile_SE.strain_detect.test.gz
echo 'running diff to compare strain_detect single-end output with expected result'
diff outfile_SE.strain_detect.gz outfile_SE.strain_detect.test.gz


echo ''
echo ''
echo 'TEST3 kmer_scrub_count'
# note we just grab the first 10K lines to keep file size small
time ../src/kmer_scrub_count -r Bacteroides_ovatus_1001283st1_B8_1001283B150210_160208.fna.gz -A scrub_strain_list.txt -B scrub_metagenome_list.txt | head -n 10000 > outfile.scrub_kmer_counts.test
echo 'running diff to compare scrub_kmer_counts output files'
diff outfile.scrub_kmer_counts outfile.scrub_kmer_counts.test
