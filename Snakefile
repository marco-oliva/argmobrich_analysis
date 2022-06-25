############################################################
## Imports
############################################################
import os

############################################################
## Messages
############################################################

############################################################
## Config file and shorthands
############################################################

configfile: "config.json"
workdir: config["WORKFLOW"]["WORKDIR"]

samples_dir = config["WORKFLOW"]["DATA"]
databases_dir = "databases"
tmp_dir = "tmp"

############################################################
## Pipeline Steps
############################################################

############################################################
# Deduplication

rule deduplicate_reads:
    input:
        reads = os.path.join(samples_dir, "{sample_name}.fastq")

    params:
        num_of_clusters = config["MISC"]["DEDUP_CLUSTERS"],
        tmp_dir_clusters = tmp_dir,
        find_duplicates_script = config["SCRIPTS"]["FIND_DUPLICATES"],
        deduplicate_script = config["SCRIPTS"]["DEDUPLICATE"],
        tmp_dir = tmp_dir

    conda:
        "envs/deduplication.yaml"
    envmodules:
        "python/3.8",
        "blat/20140318"

    threads: config["MISC"]["DEDUP_THREADS"]

    output:
        duplicates_csv = os.path.join(tmp_dir, "{sample_name}.fastq" + config["EXTENSION"]["DUPLICATES"]),
        sample_name = "{sample_name}.fastq" + config["EXTENSION"]["DEDUPLICATED"]

    shell:
        """
        mkdir -p {params.tmp_dir}
        python3 {params.find_duplicates_script} -r {input.reads} -o {params.tmp_dir_clusters} \
            -n {params.num_of_clusters} -t {threads} -b blat > {output.duplicates_csv}
        python3 {params.deduplicate_script} -r {input.reads} -d {output.duplicates_csv} > {output.sample_name}
        """

############################################################s
# Alignment

rule align_to_megares:
    input:
        reads = os.path.join(samples_dir, "{sample_name}.fastq"),
        megares_v2_seqs = os.path.join(databases_dir,"/megares_full_database_v2.00.fasta")

    params:
        minimap_flags = config["MINIMAP2"]["ALIGNER_PB_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_ONT_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_HIFI_OPTION"]

    conda:
        "envs/alignment.yaml"
    envmodules:
        "python/3.8",
        "minimap/2.21"

    threads: config["MINIMAP2"]["THREADS"]

    output:
        megares_out_sam = "{sample_name}.fastq" + config["EXTENSION"]["A_TO_MEGARES"]

    shell:
        """
        minimap2 -t {threads} {params.minimap_flags} {input.megares_v2_seqs} {input.reads} -o {output.megares_out_sam}
        """

rule align_to_mges:
    input:
        reads = os.path.join(samples_dir, "{sample_name}.fastq"),
        mges_database = os.path.join(databases_dir,"mges_combined.fa")
    params:
        minimap_flags = config["MINIMAP2"]["ALIGNER_PB_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_ONT_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_HIFI_OPTION"]

    conda:
        "envs/alignment.yaml"
    envmodules:
        "python/3.8",
        "minimap/2.21"

    threads: config["MINIMAP2"]["THREADS"]

    output:
        mges_out_sam = "{sample_name}.fastq" + config["EXTENSION"]["A_TO_MGES"]

    shell:
        """
        minimap2 -t {threads} {params.minimap_flags} {input.mges_database} {input.reads} -o {output.mges_out_sam}
        """

rule align_to_kegg:
    input:
        reads = os.path.join(samples_dir, "{sample_name}.fastq"),
        kegg_database = os.path.join(databases_dir,"kegg_genes.fa")

    params:
        minimap_flags = config["MINIMAP2"]["ALIGNER_PB_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_ONT_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_HIFI_OPTION"]

    conda:
        "envs/alignment.yaml"
    envmodules:
        "python/3.8",
        "minimap/2.21"

    threads: config["MINIMAP2"]["THREADS"]

    output:
        kegg_out_sam = "{sample_name}.fastq" + config["EXTENSION"]["A_TO_KEGG"]

    shell:
        """
        minimap2 -t {threads} {params.minimap_flags} {input.kegg_database} {input.reads} -o {output.kegg_out_sam}
        """

rule pass_config_file:
    output:
        out_config_file = "config.ini"

    run:
        import configparser
        with open(output.out_config_file,'w') as configfile_out:
            print(config.items())
            config_parser = configparser.ConfigParser()
            config_parser.read_dict(config)
            config_parser.write(configfile_out)

rule resisome_and_mobilome:
    input:
        megares_sam = "{sample_name}.fastq" + config["EXTENSION"]["A_TO_MEGARES"],
        mges_sam = "{sample_name}.fastq" + config["EXTENSION"]["A_TO_MGES"],
        config_file = "config.ini"

    params:
        resistome_mobilome_script = config["SCRIPTS"]["GEN_RESISTOME_AND_MOBILOME"],
        output_prefix = "{sample_name}.fastq_"

    conda:
        "envs/pipeline.yaml"
    envmodules:
        "python/3.8"

    output:
        resistome_richness = "{sample_name}.fastq_" + config["MISC"]["RESISTOME_STRATEGY"]
                                                    + config["EXTENSION"]["RESISTOME_RICHNESS"],
        resistome_diversity = "{sample_name}.fastq_" + config["MISC"]["RESISTOME_STRATEGY"]
                                                     + config["EXTENSION"]["RESISTOME_DIVERSITY"],
        mobilome = "{sample_name}.fastq_" + config["MISC"]["MOBILOME_STRATEGY"] + config["EXTENSION"]["MOBILOME"]

    shell:
        """
        python3 {params.resistome_mobilome_script} \
            -r {wildcards.sample_name}.fastq \
            -a {input.megares_sam} \
            -m {input.mges_sam} \
            -c {input.config_file} \
            -o {params.output_prefix}
        """

rule find_colocalizations:
    input:
        megares_sam = "{sample_name}.fastq" + config["EXTENSION"]["A_TO_MEGARES"],
        mges_sam = "{sample_name}.fastq" + config["EXTENSION"]["A_TO_MGES"],
        kegg_sam = "{sample_name}.fastq" + config["EXTENSION"]["A_TO_KEGG"],
        config_file = "config.ini"

    params:
        find_colocalizations_script = config["SCRIPTS"]["FIND_COLOCALIZATIONS"],
        output_directory = os.getcwd()

    conda:
        "envs/pipeline.yaml"
    envmodules:
        "python/3.8"

    output:
        colocalizations = "{sample_name}.fastq_" + config["EXTENSION"]["COLOCALIZATIONS"]

    shell:
        """
        python3 {params.find_colocalizations_script} \
            -r {wildcards.sample_name}.fastq \
            -a {input.megares_sam} \
            -m {input.mges_sam} \
            -k {input.kegg_sam} \
            -c {input.config_file} \
            -o {params.output_directory} \
            > {output.colocalizations}
        """

rule colocalization_richness:
    input:
        colocalizations = "{sample_name}.fastq_" + config["EXTENSION"]["COLOCALIZATIONS"],
        config_file = "config.ini"

    params:
        find_colocalizations_script = config["SCRIPTS"]["COLOCALIZATIONS_RICHNESS"]

    conda:
        "envs/pipeline.yaml"
    envmodules:
        "python/3.8"

    output:
        colocalizations_richness = "{sample_name}.fastq_" + config["EXTENSION"]["COLOCALIZATIONS"]

    shell:
        """
        python3 {params.find_colocalizations_script} \
            -i {input.colocalizations} \
            -c {input.config_file} \
            > {output.colocalizations_richness}
        """


############################################################
## Databases
############################################################

rule get_megares_v2:
    output:
        megares_v2_seqs = os.path.join(databases_dir,"megares_full_database_v2.00.fasta"),
        megares_v2_ontology = os.path.join(databases_dir,"megares_full_annotations_v2.00.csv")

    shell:
        """
        mkdir -p {databases_dir}
        wget http://megares.meglab.org/download/megares_v2.00/megares_full_database_v2.00.fasta -O {output.megares_v2_seqs}
        wget http://megares.meglab.org/download/megares_v2.00/megares_full_annotations_v2.00.csv -O {output.megares_v2_ontology}
        """

rule get_MGEs_DBs:
    output:
        mges_combined_db = os.path.join(databases_dir,"mges_combined.fa")

    params:
        mges_path = config["DATABASE"]["MGES"]

    shell:
        """
        mkdir -p {databases_dir}
        cp {params.mges_path} {output.mges_combined_db}
        """

rule get_KEGG_Prokaryotes_DBs:
    output:
        kegg_prokaryotes_db = os.path.join(databases_dir,"kegg_genes.fa")

    params:
        kegg_path = config["DATABASE"]["KEGG"]

    shell:
        """
        mkdir -p {databases_dir}
        cp {params.kegg_path} {output.kegg_prokaryotes_db}
        """

############################################################
## Cleans
############################################################

rule clean:
    shell:
        """
        rm -rf {databases_dir} {tmp_dir}
        """

rule clean_sam_files:
    shell:
        """
        rm -rf *.sam
        """
