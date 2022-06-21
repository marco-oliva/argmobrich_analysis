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

configfile: "config.yaml"

work_dir = "workdir"
databases_dir = "databases"
tmp_dir = "tmp"
tools_dir = "tools"


############################################################
## Pipeline Steps
############################################################

############################################################
# Deduplication

rule deduplicate_reads:
    input:
        reads = "{sample_name}.fastq",
        blat_sif = os.path.join(tools_dir, "blat_sif")

    params:
        num_of_clusters = config["MISC"]["DEDUP_CLUSTERS"],
        tmp_dir_clusters = tmp_dir,
        find_duplicates_script = config["SCRIPTS"]["FIND_DUPLICATES"],
        deduplicate_script = config["SCRIPTS"]["DEDUPLICATE"],
        tmp_dir = tmp_dir,
        work_dir = work_dir

    threads: config["MISC"]["DEDUP_THREADS"]

    output:
        duplicates_csv = os.path.join(tmp_dir, "{sample_name}.fastq" + config["EXTENSION"]["DUPLICATES"]),
        sample_name = os.path.join(work_dir, "{sample_name}.fastq" + config["EXTENSION"]["DEDUPLICATED"])

    shell:
        """
        mkdir -p {params.tmp_dir} {params.work_dir}
        python3 {params.find_duplicates_script} -r {input.reads} -o {params.tmp_dir_clusters} \
            -n {params.num_of_clusters} -t {threads} -b {input.blat_sif} > {output.duplicates_csv}
        python3 {params.deduplicate_script} -r {input.reads} -d {output.duplicates_csv} > {output.sample_name}
        """

############################################################s
# Alignment

rule align_to_megares:
    input:
        reads = "{sample_name}.fastq",
        megares_v2_seqs = os.path.join(databases_dir,"/megares_full_database_v2.00.fasta"),
        minimap2_sif = os.path.join(tools_dir, "minimap2_sif")

    params:
        minimap_flags = config["MINIMAP2"]["ALIGNER_PB_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_ONT_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_HIFI_OPTION"]

    output:
        megares_out_sam = os.path.join(work_dir, "{sample_name}.fastq" + config["EXTENSION"]["A_TO_MEGARES"])

    threads: config["MINIMAP2"]["THREADS"]

    shell:
        """
        {input.minimap2_sif} -t {threads} {params.minimap_flags} {input.megares_v2_seqs} {input.reads} -o {output.megares_out_sam}
        """

rule align_to_mges:
    input:
        reads = "{sample_name}.fastq",
        mges_database = os.path.join(databases_dir,"mges_combined.fa"),
        minimap2_sif = os.path.join(tools_dir, "minimap2_sif")

    params:
        minimap_flags = config["MINIMAP2"]["ALIGNER_PB_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_ONT_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_HIFI_OPTION"]

    output:
        mges_out_sam = os.path.join(work_dir, "{sample_name}.fastq" + config["EXTENSION"]["A_TO_MGES"])

    threads: config["MINIMAP2"]["THREADS"]

    shell:
        """
        {input.minimap2_sif} -t {threads} {params.minimap_flags} {input.mges_database} {input.reads} -o {output.mges_out_sam}
        """

rule align_to_kegg:
    input:
        reads = "{sample_name}.fastq",
        kegg_database = os.path.join(databases_dir,"kegg_genes.fa"),
        minimap2_sif = os.path.join(tools_dir, "minimap2_sif")

    params:
        minimap_flags = config["MINIMAP2"]["ALIGNER_PB_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_ONT_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_HIFI_OPTION"]

    output:
        kegg_out_sam = os.path.join(work_dir, "{sample_name}.fastq" + config["EXTENSION"]["A_TO_KEGG"])

    threads: config["MINIMAP2"]["THREADS"]

    shell:
        """
        {input.minimap2_sif} -t {threads} {params.minimap_flags} {input.kegg_database} {input.reads} -o {output.kegg_out_sam}
        """

rule produce_config_file:
    output:
        out_config_file = os.path.join(work_dir, "config.ini")

    run:
        import configparser
        with open(output.out_config_file,'w') as configfile_out:
            config.write(configfile_out)

rule resisome_and_mobilome:
    input:
        megares_sam = "{sample_name}.fastq" + config["EXTENSION"]["A_TO_MEGARES"],
        mges_sam = "{sample_name}.fastq" + config["EXTENSION"]["A_TO_MGES"]

    params:
        resistome_mobilome_script = config["SCRIPTS"]["GEN_RESISTOME_AND_MOBILOME"]

    output:
        resistome_richness = os.path.join(work_dir,"{sample_name}.fastq_" + config["MISC"]["RESISTOME_STRATEGY"]
                                                    + config["EXTENSION"]["RESISTOME_RICHNESS"]),
        resistome_diversity = os.path.join(work_dir,"{sample_name}.fastq_" + config["MISC"]["RESISTOME_STRATEGY"]
                                                    + config["EXTENSION"]["RESISTOME_DIVERSITY"]),
        mobilome = os.path.join(work_dir,"{sample_name}.fastq_" + config["MISC"]["MOBILOME_STRATEGY"]
                                                    + config["EXTENSION"]["MOBILOME"])

    shell:
        """
        python3 {params.resistome_mobilome_script} \
            -r {wildcards.sample_name}.fastq
        """


############################################################
## Databases
############################################################

rule get_megares_v2:
    output:
        megares_v2_seqs = os.path.join(databases_dir,"/megares_full_database_v2.00.fasta"),
        megares_v2_ontology = os.path.join(databases_dir,"megares_full_annotations_v2.00.csv")

    params:
        dbs_dir = databases_dir

    shell:
        """
        mkdir -p {params.dbs_dir}
        wget http://megares.meglab.org/download/megares_v2.00/megares_full_database_v2.00.fasta -O {output.megares_v2_seqs}
        wget http://megares.meglab.org/download/megares_v2.00/megares_full_annotations_v2.00.csv -O {output.megares_v2_ontology}
        """

rule get_MGEs_DBs:
    output:
        mges_combined_db = os.path.join(databases_dir,"mges_combined.fa")

    params:
        mges_path = config["DATABASE"]["MGES"],
        dbs_dir = databases_dir

    shell:
        """
        mkdir -p {params.dbs_dir}
        cp {params.mges_path} {output.mges_combined_db}
        """

rule get_KEGG_Prokaryotes_DBs:
    output:
        kegg_prokaryotes_db = os.path.join(databases_dir,"kegg_genes.fa")

    params:
        kegg_path = config["DATABASE"]["KEGG"],
        dbs_dir = databases_dir

    shell:
        """
        mkdir -p {params.dbs_dir}
        cp {params.kegg_path} {output.kegg_prokaryotes_db}
        """

############################################################
## Tools
############################################################

rule container_minimap2:
    output:
        sif = os.path.join(tools_dir,"minimap2_sif")

    params:
        container_link = config["CONTAINERS"]["MINIMAP2"],
        tls_dir = tools_dir

    shell:
        """
        module load singularity
        mkdir -p {params.tls_dir}
        singularity pull {output.sif} docker://{params.container_link}
        """


rule container_blat:
    output:
        sif = os.path.join(tools_dir,"blat_sif")

    params:
        container_link = config["CONTAINERS"]["BLAT"],
        tls_dir = tools_dir

    shell:
        """
        module load singularity
        mkdir -p {params.tls_dir}
        singularity pull {output.sif} docker://{params.container_link}
        """


############################################################
## Cleans
############################################################

rule clean_sam_files:
    shell:
        """
        rm -rf {work_dir}/*.sam
        """

rule clean_work_dir:
    shell:
        """
        rm -rf {work_dir}
        """

rule clean_tools:
    shell:
        """
        rm -rf {tools_dir}
        """

rule clean_all:
    shell:
        """
        rm -rf {work_dir} {databases_dir} {tmp_dir} {tools_dir}
        """
