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
```

make database 

```bash
diamond makedb --in nr.fa -d nr --taxonmap ~/projects/diamond-test/prot.accession2taxid.gz --taxonnodes ~/projects/diamond-test/nodes.dmp
```

Now we can run diamond using default parameters

```bash
ln -s ../04-synthetic_read_data_processing/samples/Slycopersicum.5k.spike.fa .
# clean up unwanted symbols
# get only tomato sequences 
cat Slycopersicum.5k.spike.fa | grep "lycopersicum" -A 1 | sed 's/-//g' > temp
mv temp Slycopersicum.5k.spike.fa
# run diamond
diamond blastx -d ~/refdb/ncbi/nr.dmnd -q Slycopersicum.5k.spike.fa -o matches.lca --outfmt 102 --taxonmap prot.accession2taxid.gz --taxonnodes nodes.dmp
```

Results (format: count taxonID)

```bash
awk '{print $2}' matches.lca | sort | uniq -c
   4847 0
      2 1
      1 107324
     30 131567
      5 1437183
      9 1437201
      1 2303987
      9 28526
      1 3398
      1 4069
     13 4070
      1 4071
      1 4072
     12 4081
     16 4083
     18 4107
      2 4113
      5 424551
     20 49274
      2 50514
      1 562
      1 58024
      1 71274
      1 91888
```
