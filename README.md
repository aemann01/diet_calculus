# diet_calculus

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4265311.svg)](https://doi.org/10.5281/zenodo.4265311)

This repository describes the analysis performed in the paper Mann _et al._ 2020. Do I have something in my teeth? The trouble with genetic analyses of diet from archaeological dental calculus (*In press*)

## Setup

This repository assumes you are running in a Unix environment (e.g. Mac OSX or
Linux) and you have conda installed.

To get this repository:

- Install and set up up anaconda or miniconda as described at the [bioconda
  documentation](https://bioconda.github.io/user/install.html), including
  setting up channels.
- Clone this repository to your machine and change into the directory with

```bash
git clone https://github.com/aemann01/diet_calculus.git && cd diet_calculus/
```

- Run the following command to install the one of the environments

```bash
conda env create -f environment.yml

```

- To load a given environment run

```bash
conda activate mann2020_<analysis_name>
```

- To turn off the environment run

```bash
conda deactivate
```

## Structure

This repository contains a directory for each analysis.

Within each directories a Conda `environment.yaml` file for the particular
analysis' software, a README describing the procedure of analysis, and
intermediate and results files. Note that some intermediate files are not
contained here due to large sizes (e.g. raw sequencing data and databases).

### 01-example_mapdamage

### 02-damage_downscaling

### 03-synthetic_dataset_generation

### 04-synthetic_dataset_processing

### 05-taxonomic_profiling

