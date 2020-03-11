#Data processing for dietary analysis of calculus paper
#Author: Allison E. Mann 2020

###################
#ENVIRONMENT SETUP
###################
mkdir samples sim mock_euks mock_oral kraken examples sim/bact sim/endo sim/cont adapterremoval

###############
#DOWNLOAD DATA
###############
echo "Downloading data"
cd mock_euks
wget -i ../euk.ids -q &
cd ../examples 
wget -i ../examples.ids -q & 
cd ../mock_oral
wget -i ../bac.ids -q
echo "Download complete"

#############################
#REMOVE CONTAMINATED CONTIGS
#############################
cd ../mock_euks
ls *gz | parallel 'gzip -d {}'
grep ">" *fna | awk -F":" '{print $2}' | sed 's/>//' | awk '{print $1}' > query.txt
comm -12 <(sort query.txt) <(sort ../contaminated_genomes.Steinegger2020.txt) > contamination_hits.txt
#how many contaminated contigs?
wc -l contamination_hits.txt
#how many contigs do we have per organism?
grep ">" *fna -c
# GCA_002197005.1_CerEla1.0_genomic.fna:11479
# GCA_900519105.1_iwgsc_refseqv1.0_genomic.fna:22
# GCF_000001405.39_GRCh38.p13_genomic.fna:639
# GCF_000002315.6_GRCg6a_genomic.fna:464
# GCF_000003025.6_Sscrofa11.1_genomic.fna:613
# GCF_000003625.3_OryCun2.0_genomic.fna:3241
# GCF_000005005.2_B73_RefGen_v4_genomic.fna:267
# GCF_000188115.4_SL3.0_genomic.fna:3150
# GCF_000233375.1_ICSASG_v2_genomic.fna:232155
# GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fna:192
# GCF_003086295.2_arahy.Tifrunner.gnm1.KYV3_genomic.fna:385
#remove contaminated seqs
ls *fna | sed 's/.fna//' | while read line; do awk 'BEGIN{while((getline<"contamination_hits.txt")>0)l[">"$1]=1}/^>/{f=!l[$1]}f' $line.fna > $line.fix.fa; done
#how much/who lost contigs?
grep ">" *fix.fa -c
# GCA_002197005.1_CerEla1.0_genomic.fix.fa:11460
# GCA_900519105.1_iwgsc_refseqv1.0_genomic.fix.fa:19
# GCF_000001405.39_GRCh38.p13_genomic.fix.fa:639
# GCF_000002315.6_GRCg6a_genomic.fix.fa:464
# GCF_000003025.6_Sscrofa11.1_genomic.fix.fa:613
# GCF_000003625.3_OryCun2.0_genomic.fix.fa:3240
# GCF_000005005.2_B73_RefGen_v4_genomic.fix.fa:267
# GCF_000188115.4_SL3.0_genomic.fix.fa:3099
# GCF_000233375.1_ICSASG_v2_genomic.fix.fa:232111
# GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fix.fa:192
# GCF_003086295.2_arahy.Tifrunner.gnm1.KYV3_genomic.fix.fa:385
#clean up
rm *fna
rename 's/.fix.fa/.fna/' *fa
#now on bacterial contigs
cd ../mock_oral 
ls *gz | parallel 'gzip -d {}'
grep ">" *fna | awk -F":" '{print $2}' | sed 's/>//' | awk '{print $1}' > query.txt
comm -12 <(sort query.txt) <(sort ../contaminated_genomes.Steinegger2020.txt) > contamination_hits.txt
#no contaminated contigs detected

###############
#TAG SEQUENCES
###############
tagSeq()
{
    set $file
    for i in $name; do
        sed "s/^>/>${i}_/" ${1} > ${1}_fix
        shift
    done
}
name=$(awk '{print $2}' ../bac.rename)
file=$(awk '{print $1}' ../bac.rename)
tagSeq
rm *fna
cat *fix > mock_oral.fa
rm *fix
cd ../mock_euks
name=$(awk '{print $2}' ../euk.rename)
file=$(awk '{print $1}' ../euk.rename)
tagSeq
#wheat has to be processed separately (12GB file)
mv GCA_900519105.1_iwgsc_refseqv1.0_genomic.fna_fix ../sim/endo/
rm *fna
cat *fix > mock_euks.fa
rm *fix
cd ..

#######################
#SIMULATED ANCIENT DNA
#######################
echo "Running Gargammel"
echo "Generating aDNA from Eukaryotic genomes"
#first generate reads for wheat
cp mock_oral/mock_oral.fa sim/bact/
cp mock_oral/mock_oral.fa sim/cont/
mv sim/endo/GCA_900519105.1_iwgsc_refseqv1.0_genomic.fna_fix sim/endo/GCA_900519105.1_iwgsc_refseqv1.0_genomic.fna
perl /home/lymelab/gargammel/gargammel.pl -n 1000000 --comp 0,0,1 --minsize 25 --maxsize 125 -damage 0.03,0.4,0.01,0.3 --loc 4.106487474 --scale 0.358874723 -o sim/wheat_sim sim
mv sim/endo/GCA_900519105.1_iwgsc_refseqv1.0_genomic.fna* mock_euks/
cp mock_euks/mock_euks.fa sim/endo/
#now generate reads for all other euk genomes
perl /home/lymelab/gargammel/gargammel.pl -n 5005000 --comp 0,0,1 --minsize 25 --maxsize 125 -damage 0.03,0.4,0.01,0.3 --loc 4.106487474 --scale 0.358874723 -o sim/mock_euks_sim sim
echo "Generating aDNA from Bacterial/Archaeal genomes"
rm sim/endo/mock_euks.fa
cp mock_oral/mock_oral.fa sim/endo
#for oral genomes
perl /home/lymelab/gargammel/gargammel.pl -n 6000000 --comp 0,0,1 --minsize 25 --maxsize 125 -damage 0.03,0.4,0.01,0.3 --loc 4.106487474 --scale 0.358874723 -o sim/mock_oral_sim sim

###########################
#QUALITY FILTER/TRIM/MERGE
###########################
echo "Adapter removal, quality filter, merge reads"
AdapterRemoval --file1 sim/wheat_sim_s1.fq.gz --file2 sim/wheat_sim_s2.fq.gz --trimns --trimqualities --minquality 25 --gzip --collapse --basename mock_wheat --minlength 25 &
AdapterRemoval --file1 sim/mock_euks_sim_s1.fq.gz --file2 sim/mock_euks_sim_s2.fq.gz --trimns --trimqualities --minquality 25 --gzip --collapse --basename mock_euks --minlength 25 &
AdapterRemoval --file1 sim/mock_oral_sim_s1.fq.gz --file2 sim/mock_oral_sim_s2.fq.gz --trimns --trimqualities --minquality 25 --gzip --collapse --basename mock_oral --minlength 25 

#for real samples first have to determine the adapter sequence
AdapterRemoval --identify-adapters --file1 examples/ELSIDRON1L7_lTACTG_rCTCGA_R1.fastq.gz --file2 examples/ELSIDRON1L7_lTACTG_rCTCGA_R2.fastq.gz
AdapterRemoval --identify-adapters --file1 examples/JAE016.A0101_R1_humfilt.fastq.gz --file2 examples/JAE016.A0101_R2_humfilt.fastq.gz
AdapterRemoval --identify-adapters --file1 examples/H10b_calc_shotgun_R1.fastq.gz --file2 examples/H10b_calc_shotgun_R2.fastq.gz
#they all have same adapters but already started run, could do this in loop to clean up
AdapterRemoval --file1 examples/ELSIDRON1L7_lTACTG_rCTCGA_R1.fastq.gz --file2 examples/ELSIDRON1L7_lTACTG_rCTCGA_R2.fastq.gz --trimns --trimqualities --minquality 25 --gzip --collapse --basename neanderthal --minlength 25 --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT &
AdapterRemoval --file1 examples/JAE016.A0101_R1_humfilt.fastq.gz --file2 examples/JAE016.A0101_R2_humfilt.fastq.gz --trimns --trimqualities --minquality 25 --gzip --collapse --basename modHuman --minlength 25 --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
ls examples/*calc*R1* | sed 's/_R1.fastq.gz//' | sed 's/examples\///' | parallel 'AdapterRemoval --file1 examples/{}_R1.fastq.gz --file2 examples/{}_R2.fastq.gz --trimns --trimqualities --minquality 25 --gzip --collapse --basename {} --minlength 25 --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'
mv *gz *settings adapterremoval
rm -r sim

##############################
#FASTQ TO FASTA/DEREPLICATION
##############################
echo "Converting to fasta, dereplication"
#keeping collapsed, singletons, and high qualtiy forward reads
ls adapterremoval/*pair1*gz | parallel 'gzip -d {}' &
ls adapterremoval/*singleton*gz | parallel 'gzip -d {}' &
ls adapterremoval/*collapsed*gz | parallel 'gzip -d {}' 
cd adapterremoval
ls | grep -v "gz" | grep -v "settings" | parallel 'fastq_to_fasta -i {} -o {}.fna'
ls *.collapsed.fna | sed 's/.collapsed.fna//' | while read line; do cat $line.collapsed.fna $line.collapsed.truncated.fna $line.pair1.fna $line.pair1.truncated.fna $line.singletons.fna $line.singletons.truncated.fna > $line.fa; done 
rm *fna 
#Dereplicate
ls *fa | sed 's/.fa//' | parallel 'vsearch --derep_fulllength {}.fa --output ../{}.uniq.fa'
cd ..
rm -r adapterremoval
mv *calc* neanderthal.uniq.fa modHuman.uniq.fa samples/

###################
#TRIM POLY A TAILS 
###################
ls *uniq* | sed 's/.uniq.fa//' | while read line; do cutadapt -a "A{100}" -o $line.trim.fa $line.uniq.fa; done
rm *uniq*

#####################
#EUK GENOME SPIKE IN
#####################
echo "Splitting mock Eukaryotes by species"
rm -r temp
mkdir temp
echo "Mock oral community prep"
#fix headers so you can retrieve after kraken
sed 's/:.*//' mock_euks.trim.fa | awk '/>/{print $0=$0"_"(++i)}!/>/' > temp.mock
mv temp.mock mock_euks.trim.fa
sed 's/:.*//' mock_oral.trim.fa | awk '/>/{print $0=$0"_"(++i)}!/>/' > temp.mock
mv temp.mock mock_oral.trim.fa
sed 's/:.*//' mock_wheat.trim.fa | awk '/>/{print $0=$0"_"(++i)}!/>/' > temp.mock
mv temp.mock mock_wheat.trim.fa

#split out by species
cat mock_wheat.trim.fa mock_euks.trim.fa > temp.seqs
mv temp.seqs mock_euks.trim.fa
awk '{print $2}' euk.rename | parallel 'grep {} -A 1 mock_euks.trim.fa > temp/{}.fa'
echo "How many dietary reads post quality filter and dereplication?"
grep ">" temp/*fa -c
#generate mock oral community files
#5k depth
cat mock_oral.trim.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 4995000 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_5k.fa
#500 depth
cat mock_oral.trim.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 4999500 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_500.fa
#50 depth
cat mock_oral.trim.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 4999950 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_50.fa
echo "Mock oral files generated. Should be: 4995000, 499500, 4999950"
grep ">" *oral*5*.fa -c
echo "Mock Eukaryotic community prep"
cd temp
ls *fa | sed 's/.fa$//' > diet.ids
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
grep ">" samples/*spike*fa -c
rm -r temp
mv neanderthal.trim.fa *calc*trim* modHuman.trim.fa samples/

## TO DO--> put clean up step here, only keep files you need

#######################
#KRAKEN DATABASE BUILD
#######################
#only run this once! takes a long time to build
# mkdir ~/reference_databases/kraken_nt 
# kraken2-build --download-taxonomy --db kraken_nt
# kraken2-build --download-library nt --db kraken_nt
# kraken2-build --build --db kraken_nt --kmer-len 45 --threads 10

########
#KRAKEN
########
cd samples 
ls *fa | parallel 'gzip {}'
ls *fa.gz | sed  's/.fa.gz//' | while read line; do kraken2 --db ~/reference_databases/kraken_nt/ --threads 8 --use-names --gzip-compressed --output ../kraken/$line.out $line.fa.gz; done

#get summary of dietary hits for each spike level output
cd ../kraken
ls *5k* | sed 's/.5k.spike.out//' | while read line; do grep $line $line.5k.spike.out | awk -F"\t" '{print $2, "\t", $3}' | sed 's/^M_//' | sed 's/^MT_//' | sed 's/_.*\t/\t/' | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/' | sed 's/ //' | sed 's/ /_/g' | sed 's/_(taxid_/\t/' | sed 's/)$//' >> 5k.spike.summary; done
ls *5k* | sed 's/.5k.spike.out//' | while read line; do grep $line $line.500.spike.out | awk -F"\t" '{print $2, "\t", $3}' | sed 's/^M_//' | sed 's/^MT_//' | sed 's/_.*\t/\t/' | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/' | sed 's/ //' | sed 's/ /_/g' | sed 's/_(taxid_/\t/' | sed 's/)$//' >> 500.spike.summary; done
ls *5k* | sed 's/.5k.spike.out//' | while read line; do grep $line $line.50.spike.out | awk -F"\t" '{print $2, "\t", $3}' | sed 's/^M_//' | sed 's/^MT_//' | sed 's/_.*\t/\t/' | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/' | sed 's/ //' | sed 's/ /_/g' | sed 's/_(taxid_/\t/' | sed 's/)$//' >> 50.spike.summary; done

#bacterial summary
awk -F"\t" '{print $2, "\t", $3}' Taestivum.50.spike.out | sed 's/^M_//' | sed 's/^MT_//' | sed 's/_.*\t/\t/' | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/' | sed 's/ //' | sed 's/ /_/g' | sed 's/_(taxid_/\t/' | sed 's/)$//' > bacterial.summary
#summary for real samples
awk -F"\t" '{print $3}' modHuman.trim.out | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/' | sed 's/ /_/g' | sed 's/_(taxid_/\t/' | sed 's/)//g' > modHuman.summary &
awk -F"\t" '{print $3}' neanderthal.trim.out | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/' | sed 's/ /_/g' | sed 's/_(taxid_/\t/' | sed 's/)//g' > neanderthal.summary &
ls *calc*trim* | sed 's/.out//' | while read line; do awk -F"\t" '{print $3}' $line.out | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/' | sed 's/ /_/g' | sed 's/_(taxid_/\t/' | sed 's/)//g' > $line.summary; done
-cd ..
ls kraken/*out | parallel 'gzip {}' & 
#add header line
ls kraken/5* | while read line; do sed -i '1 i\count\tquery\tsubject\ttaxid' $line; done &
sed -i '1 i\count\tquery\tsubject\ttaxid' kraken/bacterial.summary
sed -i -e '1 i\count\tsubject\ttaxid' -e 's/^[ \t]*//' -e 's/ /\t/' kraken/modHuman.summary
sed -i -e '1 i\count\tsubject\ttaxid' -e 's/^[ \t]*//' -e 's/ /\t/' kraken/neanderthal.summary
ls kraken/*calc*summary | while read line; do sed -i -e '1 i\count\tsubject\ttaxid' -e 's/^[ \t]*//' -e 's/ /\t/' $line; done

#merge with taxonomy string
cd kraken
gzip -d ../taxid_taxonomystr.txt.gz 
ls *summary | sed 's/.summary//' | parallel 'python ../merge_tax.py -i {}.summary -o {}.merged -t ../taxid_taxonomystr.txt'
# rm *summary

