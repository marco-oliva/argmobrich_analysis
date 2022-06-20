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

############################################################
# Alignment

rule align_to_megares:
    input:
        sample_name = "{sample_name}.fastq",
        megares_v2_seqs = os.path.join(databases_dir,"/megares_full_database_v2.00.fasta"),
        minimap2_sif = os.path.join(tools_dir, "minimap2_sif")

    params:
        minimap_flags = config["MINIMAP2"]["ALIGNER_PB_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_ONT_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_HIFI_OPTION"]

    output:
        megares_out_sam = os.path.join(work_dir, "{sample_name}" + config["EXTENSION"]["A_TO_MEGARES"])

    threads: config["MINIMAP2"]["THREADS"]

    shell:
        """
        {input.minimap2_sif} -t {threads} {params.minimap_flags} {input.megares_v2_seqs} {input.sample_name} -o {output.megares_out_sam}
        """

rule align_to_mges:
    input:
        sample_name = "{sample_name}.fastq",
        mges_database = os.path.join(databases_dir,"mges_combined.fa"),
        minimap2_sif = os.path.join(tools_dir, "minimap2_sif")

    params:
        minimap_flags = config["MINIMAP2"]["ALIGNER_PB_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_ONT_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_HIFI_OPTION"]

    output:
        mges_out_sam = os.path.join(work_dir, "{sample_name}" + config["EXTENSION"]["A_TO_MGES"])

    threads: config["MINIMAP2"]["THREADS"]

    shell:
        """
        {input.minimap2_sif} -t {threads} {params.minimap_flags} {input.mges_database} {input.sample_name} -o {output.mges_out_sam}
        """

rule align_to_kegg:
    input:
        sample_name = "{sample_name}.fastq",
        kegg_database = os.path.join(databases_dir,"kegg_genes.fa"),
        minimap2_sif = os.path.join(tools_dir, "minimap2_sif")

    params:
        minimap_flags = config["MINIMAP2"]["ALIGNER_PB_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_ONT_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_HIFI_OPTION"]

    output:
        kegg_out_sam = os.path.join(work_dir, "{sample_name}" + config["EXTENSION"]["A_TO_KEGG"])

    threads: config["MINIMAP2"]["THREADS"]

    shell:
        """
        {input.minimap2_sif} -t {threads} {params.minimap_flags} {input.kegg_database} {input.sample_name} -o {output.kegg_out_sam}
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

