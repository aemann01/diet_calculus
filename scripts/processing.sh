###########################
#SIMULATED DATA PROCESSING
###########################
cd /home/lymelab/lab_members/mann/
#generate simulated bacterial, human, and dietary spikeins
perl ~/gargammel/gargammel.pl -n 50000000 --comp 0,0,1 --minsize 25 --maxsize 125 -damage 0.03,0.4,0.01,0.3 --loc 4.106487474 --scale 0.358874723 -o diet_calculus/mock_oral_simulated diet_calculus
#concatenate all eukaryotes and run same command, split out by species for individual tests
cat bact/GC*fna cont/GCF_000001405.39_GRCh38.p13_genomic_human.fna > all_euks.fna
mv all_euks.fna endo
mv endo/mock_oral.fna* .
perl ~/gargammel/gargammel.pl -n 50000000 --comp 0,0,1 --minsize 25 --maxsize 125 -damage 0.03,0.4,0.01,0.3 --loc 4.106487474 --scale 0.358874723 -o diet_calculus/mock_euks_simulated diet_calculus
