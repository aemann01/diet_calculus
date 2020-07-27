# Do I have something in my teeth? The trouble with genetic analyses of diet from archaeological dental calculus

Allison E. Mann, James A. Fellows Yates, Zandra Fagernäs, Rita M. Austin, Elizabeth A. Nelson, Courtney A. Hofman

This repository describes the analysis performed in the paper Mann _et al._ XXX.

## Setup

We provide a [conda](https://docs.conda.io/en/latest/) environment to allow reproduction of the analysis recorded in this repository, by ensuring the same software versions are used when running the commands as described in `scripts/`.

To set up this environment:

- Install and set up up anaconda or miniconda as described at the [bioconda documentation](https://bioconda.github.io/user/install.html), including setting up channels.
- Clone this repository to your machine and change into the directory with

```bash
git clone https://github.com/aemann01/diet_calculus.git && cd diet_calculus/
```

- Run the following command to install the environment

```bash
conda env create -f environment.yml

```

- To load the environment run 

```bash
conda activate qi_diet_dna_calculus
```
