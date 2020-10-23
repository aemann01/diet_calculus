To demonstrate that the patterns of misidentification seen in our kraken run are not pipeline specific we also tested the best performing diet species, tomato, using diamond (blastx).

Install diamond using conda

```bash
 conda install -c bioconda diamond 
 ```

First need to download and format the nr protein database from ncbi

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
gzip -d nr.gz
mv nr nr.fa
```

Get taxonomy files

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
unzip taxdmp.zip &
gzip -d prot.accession2taxid.gz
```

make database 

```bash
diamond makedb --in nr.fa -d nr --taxonmap ~/projects/diamond-test/prot.accession2taxid --taxonnodes ~/projects/diamond-test/nodes.dmp
```

Now we can run diamond using default parameters

```bash
ln -s ../04-synthetic_read_data_processing/samples/Slycopersicum.5k.spike.fa .
# clean up unwanted symbols
# get only tomato sequences 
cat Slycopersicum.5k.spike.fa | grep "lycopersicum" -A 1 | sed 's/-//g' > temp
mv temp Slycopersicum.5k.spike.fa
# run diamond
diamond blastx -d nr -q Slycopersicum.5k.spike.fa -o matches.lca --outfmt 102
```

Next, does retrieval of reads improve when using a more curated database? Here we will run kraken using a smaller version of NCBI's plant RefSeq database 

Make kraken database

```bash
kraken2-build --download-taxonomy --db kraken_plant
kraken2-build --download-library plant --db kraken_plant
kraken2-build --build --db kraken_plant --kmer-len 45 --threads 10
```

Zip up fasta file

```bash
gzip Slycopersicum.5k.spike.fa.gz
```

Process tomato (5k) sample

```bash
kraken2 --db kraken_plasmid/ --threads 8 --use-names --gzip-compressed --output Slycopersicum.out temp.gz
```

Results

```bash
