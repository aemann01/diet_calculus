#Data processing for dietary analysis of calculus paper

cd /home/lymelab/lab_members/mann/diet_calculus

#######################
#SIMULATED ANCIENT DNA
#######################
#generate mock oral community
perl /home/lymelab/gargammel/gargammel.pl -n 6000000 --comp 0,0,1 --minsize 25 --maxsize 125 -damage 0.03,0.4,0.01,0.3 --loc 4.106487474 --scale 0.358874723 -o sim/mock_oral_simulated sim
#generate mock dietary reads
mv mock_oral/mock_oral_raw* sim/endo/
mv sim/endo/mock_oral_raw.fa* mock_oral/
mv all_euks.fna* sim/endo/
perl /home/lymelab/gargammel/gargammel.pl -n 6000000 --comp 0,0,1 --minsize 25 --maxsize 125 -damage 0.03,0.4,0.01,0.3 --loc 4.106487474 --scale 0.358874723 -o sim/mock_euks_simulated sim

######################
#QUALITY FILTER/MERGE
######################
AdapterRemoval --file1 sim/mock_oral_simulated_s1.fq.gz --file2 sim/mock_oral_simulated_s2.fq.gz --trimns --trimqualities --minquality 25 --gzip --collapse --basename mock_oral --minlength 25
AdapterRemoval --file1 sim/mock_euks_simulated_s1.fq.gz --file2 sim/mock_euks_simulated_s2.fq.gz --trimns --trimqualities --minquality 25 --gzip --collapse --basename mock_euks --minlength 25 &
#weyrich data
AdapterRemoval --file1 weyrich2017/ELSIDRON1L7_lTACTG_rCTCGA_R1.fastq --file2 weyrich2017/ELSIDRON1L7_lTACTG_rCTCGA_R2.fastq --trimns --trimqualities --minquality 25 --gzip --collapse --basename elsidron1l7 --minlength 25 &
#velsko data
AdapterRemoval --file1 velsko2019/JAE016.A0101_R1_humfilt.fastq.gz --file2 velsko2019/JAE016.A0101_R2_humfilt.fastq.gz --trimns --trimqualities --minquality 25 --gzip --collapse --basename elsidron1l7 --minlength 25 &
rm -r adapterremoval
mkdir adapterremoval
mv *gz *settings adapterremoval

#####################
#EUK GENOME SPIKE IN
#####################
#fastq to fasta -- keeping only collapsed and collapsed truncated reads
ls adapterremoval/*collapsed* | parallel 'gzip -d {}'
ls adapterremoval/*collapsed* | parallel 'fastq_to_fasta -i {} -o {}.fa'
#split euk file by species
cat adapterremoval/mock_euks.collapsed.fa adapterremoval/mock_euks.collapsed.truncated.fa > mock_euks.fa
cat adapterremoval/mock_oral.collapsed.fa adapterremoval/mock_oral.collapsed.truncated.fa > mock_oral.fa
grep ">" mock_euks.fa | awk -F"_" '{print $2}' | sort | uniq | parallel 'grep {} -A 1 mock_euks.fa > diet_reads/{}.fa'
grep ">" diet_reads/*fa -c
# chicken.fa:156679
# corn.fa:312814
# human.fa:462216
# mock_euks.fa:4999704
# mock_oral.fa:4999747
# peach.fa:33439
# peanut.fa:379257
# pig.fa:367057
# rabbit.fa:386155
# red.fa:287740
# salmon.fa:386693
# tomato.fa:111182
# wheat.fa:2116472

##########################
#MOCK ORAL COMMUNITY PREP
##########################
#5000 depth
cat mock_oral.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 4995000 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_5k.fa
#500 depth
cat mock_oral.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 4999500 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_500.fa
#50 depth
cat mock_oral.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 4999950 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_50.fa

##############
#MOCK SAMPLES
##############
cd diet_reads
ls *fa > diet.ids
#5000 subsample
cat diet.ids | sed 's/.fa//' | while read line; do cat $line.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 5000 | awk '{printf("%s\n%s\n",$1,$2)}' > $line.5k.fa; done
#500 subsample
cat diet.ids | sed 's/.fa//' | while read line; do cat $line.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 500 | awk '{printf("%s\n%s\n",$1,$2)}' > $line.500.fa; done
#50 subsample
cat diet.ids | sed 's/.fa//' | while read line; do cat $line.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 50 | awk '{printf("%s\n%s\n",$1,$2)}' > $line.50.fa; done
#concatenate appropriate sets, assign sample names
cat diet.ids | sed 's/.fa//' | while read line; do cat $line.5k.fa ../mock_oral_5k.fa > ../sim_samples/$line.5k.spike.fa; done
cat diet.ids | sed 's/.fa//' | while read line; do cat $line.500.fa ../mock_oral_500.fa > ../sim_samples/$line.500.spike.fa; done
cat diet.ids | sed 's/.fa//' | while read line; do cat $line.50.fa ../mock_oral_50.fa > ../sim_samples/$line.50.spike.fa; done
cd ..
grep ">" sim_samples/*fa -c
#should all be 5 million




#######################
#KRAKEN DATABASE BUILD
#######################
#only run this once! takes a long time to build
kraken2-build --download-taxonomy --db /home/lymelab/reference_databases/kraken_nt
kraken2-build --download-library nt --db /home/lymelab/reference_databases/kraken_nt
kraken2-build --build --db /home/lymelab/reference_databases/kraken_nt --kmer-len 45 --threads 10

########
#KRAKEN
########



























