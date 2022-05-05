# Target-enriched long-read sequencing (TELS) contextualizes antibiotic resistance in microbiomes 

This repository contains the scripts used to perform the analysis descrived in the *Bioinformatics Analysis* section.


## Install

```bash
git clone https://github.com/marco-oliva/argmobrich_analysis.git
cd argmobrich_analysis
pip install -r requirements.txt
```

## Usage exmplae

Edit the configuration file `confing.ini` specifying where to find the database and the aligner to use. It is also possible to specify what steps of the pipeline to execute. Once all the details have been specified in `confing.ini` you can run the pipeline with: 

```bash
./tels.py -c config.ini -o out_dir -t 16 -i reads.fastq 
```