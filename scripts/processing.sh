#Data processing for dietary analysis of calculus paper

###################
#ENVIRONMENT SETUP
###################
mkdir samples sim mock_euks mock_oral kraken1 examples sim/bact sim/endo sim/cont adapterremoval

###############
#DOWNLOAD DATA
###############
echo "Downloading data"
cd mock_euks
wget -i ../euk.ids -q &
cd ../examples 
wget -i ../examples.ids -q & ## TO DO--> direct download weyrich data
cd ../mock_oral
wget -i ../bac.ids -q
echo "Download complete"

###############
#TAG SEQUENCES
###############
name=$(awk '{print $2}' test.txt)
file=$(awk '{print $1}' test.txt)
tagSeq()
{
    set $file
    for i in $name; do
        sed -i "s/^>.*/>${i}/" ${1}
        shift
    done
}
tagSeq
cat *fna > mock_oral.fa
rm GCF_00*fna
#now run on eukaryotic genomes
cd ../mock_euks
name=$(awk '{print $2}' test.txt)
file=$(awk '{print $1}' test.txt)
tagSeq
ls *gz | parallel 'gzip -d {}'
# wheat has to be processed separately -- 14G genome compared to 1-2G
mv GCA_900519105.1_iwgsc_refseqv1.0_genomic.fna ../sim/endo/
#everyone else concatenated into single fasta
cat *fna > mock_euks.fa
rm GC*fna &
cd ..

#######################
#SIMULATED ANCIENT DNA
#######################
echo "Running Gargammel"
echo "Generating aDNA from Eukaryotic genomes"
#first generate reads for wheat
cp mock_oral/mock_oral.fa sim/bact/
cp mock_oral/mock_oral.fa sim/cont/
perl /home/lymelab/gargammel/gargammel.pl -n 1000000 --comp 0,0,1 --minsize 25 --maxsize 125 -damage 0.03,0.4,0.01,0.3 --loc 4.106487474 --scale 0.358874723 -o sim/wheat_sim sim
rm /sim/endo/GCA_900519105.1_iwgsc_refseqv1.0_genomic.fna*
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
AdapterRemoval --file1 examples/ELSIDRON1L7_lTACTG_rCTCGA_R1.fastq.gz --file2 examples/ELSIDRON1L7_lTACTG_rCTCGA_R2.fastq.gz --trimns --trimqualities --minquality 25 --gzip --collapse --basename neanderthal --minlength 25 &
AdapterRemoval --file1 examples/JAE016.A0101_R1_humfilt.fastq.gz --file2 examples/JAE016.A0101_R2_humfilt.fastq.gz --trimns --trimqualities --minquality 25 --gzip --collapse --basename modHuman --minlength 25 
AdapterRemoval --file1 sim/mock_oral_sim_s1.fq.gz --file2 sim/mock_oral_sim_s2.fq.gz --trimns --trimqualities --minquality 25 --gzip --collapse --basename mock_oral --minlength 25 
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
ls *fa | sed 's/.fa//' | parallel 'vsearch --derep_fulllength {}.fa --output ../{}.uniq.fa'
cd ..
rm -r adapterremoval
mv neanderthal.uniq.fa modHuman.uniq.fa wheat_sim.uniq.fa samples

#####################
#EUK GENOME SPIKE IN
#####################
echo "Splitting mock Eukaryotes by species"
rm -r temp
mkdir temp
grep ">" mock_euks.uniq.fa | awk -F":" '{print $1}' | sed 's/>.*_//' | sed 's/>//' | sort | uniq | parallel 'grep {} -A 1 mock_euks.uniq.fa > temp/{}.fa'
echo "How many dietary reads post quality filter and dereplication?"
grep ">" temp/*fa -c
echo "Mock oral community prep"
#5k depth
cat mock_oral.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 4995000 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_5k.fa
#500 depth
cat mock_oral.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 4999500 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_500.fa
#50 depth
cat mock_oral.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 4999950 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_50.fa
echo "Mock oral files generated. Should be: 4995000, 499500, 4999950"
grep ">" *oral*5*.fa -c
echo "Mock Eukaryotic community prep"
mv mock_wheat.uniq.fa temp
mv temp/mock_wheat.uniq.fa temp/Taestivum.fa
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

## TO DO--> put clean up step here, only keep files you need

#######################
#KRAKEN DATABASE BUILD
#######################
#only run this once! takes a long time to build
# mkdir kraken_ng 
# kraken2-build --download-taxonomy --db kraken_nt
# kraken2-build --download-library nt --db kraken_nt
# kraken2-build --build --db kraken_nt --kmer-len 45 --threads 10

########
#KRAKEN
########
cd samples 
ls *fa | parallel 'gzip {}'
ls *fa.gz | sed  's/.fa.gz//' | while read line; do kraken2 --db ~/reference_databases/kraken_nt/ --threads 8 --use-names --gzip-compressed --output ../kraken1/$line.out $line.fa.gz; done
#run 2
mkdir ../kraken2
ls *spike*fa.gz | sed  's/.fa.gz//' | while read line; do kraken2 --db ~/reference_databases/kraken_nt/ --threads 8 --use-names --gzip-compressed --output ../kraken2/$line.out $line.fa.gz; done
#run 3
mkdir ../kraken3 
ls *spike*fa.gz | sed  's/.fa.gz//' | while read line; do kraken2 --db ~/reference_databases/kraken_nt/ --threads 8 --use-names --gzip-compressed --output ../kraken3/$line.out $line.fa.gz; done

#get summaries and remove kraken output
cd ../kraken1
ls *spike*out | sed 's/.spike.out//' | while read line; do awk -F"\t" '{print $2, "\t", $3}' $line.spike.out | sed 's/:.:.*\t/\t/' | sed 's/^M_//' | sed 's/^MT_//' | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/' | sed 's/ //' | sed 's/ /_/g' | sed 's/_(taxid_/\t/' | sed 's/)$//' > $line.summary; done
#just want to get taxonomic information from neanderthal and modern human
awk -F"\t" '{print $3}' modHuman.out | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/' | sed 's/ /_/g' > modHuman.summary
awk -F"\t" '{print $3}' neanderthal.out | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/' | sed 's/ /_/g' > neanderthal.summary
cd ../kraken2
ls *out | sed 's/.spike.out//' | while read line; do awk -F"\t" '{print $2, "\t", $3}' $line.spike.out | sed 's/:.:.*\t/\t/' | sed 's/^M_//' | sed 's/^MT_//' | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/' | sed 's/ //' | sed 's/ /_/g' | sed 's/_(taxid_/\t/' | sed 's/)$//' > $line.summary; done
cd ../kraken3
ls *out | sed 's/.spike.out//' | while read line; do awk -F"\t" '{print $2, "\t", $3}' $line.spike.out | sed 's/:.:.*\t/\t/' | sed 's/^M_//' | sed 's/^MT_//' | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/' | sed 's/ //' | sed 's/ /_/g' | sed 's/_(taxid_/\t/' | sed 's/)$//' > $line.summary; done
cd ..
ls kraken1/*out | parallel 'gzip {}' & 
ls kraken2/*out | parallel 'gzip {}' & 
ls kraken3/*out | parallel 'gzip {}' & 

#add header line
ls kraken*/*5* | while read line; do sed -i '1 i\count\tquery\tsubject\ttaxid' $line; done
sed -i -e '1 i\count\tsubject\ttaxid' -e 's/^[ \t]*//' -e 's/ /\t/' kraken1/modHuman.summary
sed -i -e '1 i\count\tsubject\ttaxid' -e 's/^[ \t]*//' -e 's/ /\t/' kraken1/neanderthal.summary
#no variation between kraken runs, continue with just the first run
rm -r kraken2 kraken3
#merge with taxonomy string
cd kraken1
ls | sed 's/.summary//' | parallel 'python ../scripts/merge_tax.py -i {}.summary -o {}.merged -t ../taxid_taxonomystr.txt'
rm *summary