# Target-enriched long-read sequencing (TELSeq) contextualizes antibiotic resistance in microbiomes.

This repository contains the scripts used to perform the analysis descrived in the *Bioinformatics Analysis* section. 

### Snakemake pipeline: [`telseq`](https://github.com/marco-oliva/telseq)

A Snakemake pipeline based on the scripts contained in this repostiory is available here: [https://github.com/marco-oliva/telseq](https://github.com/marco-oliva/telseq).

[`telseq`](https://github.com/marco-oliva/telseq) is currently better supported and will be extended with additional features.

---

### Requirements

In order to run the scripts the following requirements are needed:

- python3
- pip

### Install

```bash
git clone https://github.com/marco-oliva/argmobrich_analysis.git
cd argmobrich_analysis
pip install -r requirements.txt
```

### Usage

```
usage: tels.py [-h] [-i INPUT_PATH] -o OUTPUT_DIR_PATH [-c CONFIG_PATH] [-t HELPER_THREADS]

Colocalizations Pipeline

optional arguments:
  -h, --help          show this help message and exit
  -i INPUT_PATH       Input file
  -o OUTPUT_DIR_PATH  Output directory
  -c CONFIG_PATH      Config file
  -t HELPER_THREADS   Helper threads

```

### Usage example

Edit the configuration file `confing.ini` specifying where to find the database and the aligner to use. It is also possible to specify what steps of the pipeline to execute. Once all the details have been specified in `confing.ini` you can run the pipeline with: 

```bash
./tels.py -c config.ini -o out_dir -t 16 -i reads.fastq 
```
