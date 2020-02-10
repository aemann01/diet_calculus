#Data processing for dietary analysis of calculus paper

cd /home/lymelab/lab_members/mann/diet_calculus
###################
#ENVIRONMENT SETUP
###################
#TO DO -- generate folders

###############
#DOWNLOAD DATA
###############
#TO DO -- wget pull data

#######################
#SIMULATED ANCIENT DNA
#######################
echo "Running Gargammel"
#generate mock oral community
perl /home/lymelab/gargammel/gargammel.pl -n 5002000 --comp 0,0,1 --minsize 25 --maxsize 125 -damage 0.03,0.4,0.01,0.3 --loc 4.106487474 --scale 0.358874723 -o sim/mock_oral_simulated sim
#generate mock dietary reads
mv mock_oral/mock_oral_raw* sim/endo/
mv sim/endo/mock_oral_raw.fa* mock_oral/
mv all_euks.fna* sim/endo/
perl /home/lymelab/gargammel/gargammel.pl -n 5002000 --comp 0,0,1 --minsize 25 --maxsize 125 -damage 0.03,0.4,0.01,0.3 --loc 4.106487474 --scale 0.358874723 -o sim/mock_euks_simulated sim
mv sim/end/all_euks* mock_euks

######################
#QUALITY FILTER/MERGE
######################
echo "Adapter removal and quality trim"
AdapterRemoval --file1 sim/mock_oral_simulated_s1.fq.gz --file2 sim/mock_oral_simulated_s2.fq.gz --trimns --trimqualities --minquality 25 --gzip --collapse --basename mock_oral --minlength 35 &
AdapterRemoval --file1 sim/mock_euks_simulated_s1.fq.gz --file2 sim/mock_euks_simulated_s2.fq.gz --trimns --trimqualities --minquality 25 --gzip --collapse --basename mock_euks --minlength 35 &
#weyrich data
AdapterRemoval --file1 weyrich2017/ELSIDRON1L7_lTACTG_rCTCGA_R1.fastq --file2 weyrich2017/ELSIDRON1L7_lTACTG_rCTCGA_R2.fastq --trimns --trimqualities --minquality 25 --gzip --collapse --basename elsidron1l7 --minlength 25 &
#velsko data
AdapterRemoval --file1 velsko2019/JAE016.A0101_R1_humfilt.fastq.gz --file2 velsko2019/JAE016.A0101_R2_humfilt.fastq.gz --trimns --trimqualities --minquality 25 --gzip --collapse --basename jae016 --minlength 25 
rm -r adapterremoval
mkdir adapterremoval
mv *gz *settings adapterremoval

#####################
#EUK GENOME SPIKE IN
#####################
#fastq to fasta -- keeping only collapsed and collapsed truncated reads
ls adapterremoval/*collapsed*gz | parallel 'gzip -d {}'
ls adapterremoval/*collapsed | parallel 'fastq_to_fasta -i {} -o {}.fa'




#rezip fastq
ls adapterremoval/*collapsed | grep -v ".fa" | parallel 'gzip {}' &
#concatenate
ls adapterremoval/*fa | sed 's/.collapsed.*//' | sort | uniq | while read line; do cat $line.collapsed.fa $line.collapsed.truncated.fa > $line.fa; done
#remove unecessary fa files
rm adapterremoval/*collapsed*fa 
#dereplicate
ls adapterremoval/*fa | sed 's/.fa//' | parallel 'vsearch --derep_fulllength {}.fa --output {}.uniq.fa'
mv adapterremoval/*uniq.fa .
rm adapterremoval/*fa

#split euk file by species
rm -r diet_reads/
mkdir diet_reads
grep ">" mock_euks.uniq.fa | awk -F"_" '{print $2}' | sort | uniq | parallel 'grep {} -A 1 mock_euks.uniq.fa > diet_reads/{}.fa'
echo "How many dietary reads post quality filter and dereplication?"
grep ">" diet_reads/*fa -c

##########################
#MOCK ORAL COMMUNITY PREP
##########################
#5000 depth
cat mock_oral.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 4995000 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_5k.fa
#500 depth
cat mock_oral.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 4999500 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_500.fa
#50 depth
cat mock_oral.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 4999950 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_50.fa
echo "Mock oral community preparation, 4995000, 499500, 4999950"
grep ">" *oral*5*.fa -c

##############
#MOCK SAMPLES
##############
mkdir samples 
mv jae016.uniq.fa elsidron1l7.uniq.fa samples 
cd diet_reads
ls *fa | sed 's/.fa//' > diet.ids
#5000 subsample
cat diet.ids | while read line; do cat $line.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 5000 | awk '{printf("%s\n%s\n",$1,$2)}' > $line.5k.fa; done &
#500 subsample
cat diet.ids | while read line; do cat $line.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 500 | awk '{printf("%s\n%s\n",$1,$2)}' > $line.500.fa; done &
#50 subsample
cat diet.ids | while read line; do cat $line.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 50 | awk '{printf("%s\n%s\n",$1,$2)}' > $line.50.fa; done
#concatenate appropriate sets, assign sample names
cat diet.ids | while read line; do cat $line.5k.fa ../mock_oral_5k.fa > ../samples/$line.5k.spike.fa; done &
cat diet.ids | while read line; do cat $line.500.fa ../mock_oral_500.fa > ../samples/$line.500.spike.fa; done &
cat diet.ids | while read line; do cat $line.50.fa ../mock_oral_50.fa > ../samples/$line.50.spike.fa; done
cd ..
echo "Generate mock samples, should all be 5000000"
grep ">" samples/*fa -c

#clean up
rm -r *5*fa *uniq*fa diet_reads adapterremoval 
ls samples/*fa | parallel 'gzip {}'

#######################
#KRAKEN DATABASE BUILD
#######################
#only run this once! takes a long time to build
# kraken2-build --download-taxonomy --db /home/lymelab/reference_databases/kraken_nt
# kraken2-build --download-library nt --db /home/lymelab/reference_databases/kraken_nt
# kraken2-build --build --db /home/lymelab/reference_databases/kraken_nt --kmer-len 45 --threads 10

########
#KRAKEN
########
mkdir kraken
cd samples 
ls *fa.gz | sed  's/.spike.fa.gz//' | while read line; do kraken2 --db ~/reference_databases/kraken_nt/ --threads 8 --unclassified-out ../kraken/$line.unclassified.out --use-names --gzip-compressed --output ../kraken/$line.out $line.spike.fa; done
























