# snakemake --cores 10 --dry-run
# snakemake --cores 10 --until wget_uniprot_fasta --dry-run

import pandas as pd

#sequence_ids = pd.read_csv("runs/unstructured_dms_tabfix.tsv", sep='\t')['Uniprot ID']

rule all:
    input:
        "results/FIKKs/FIKK3/model/model_1.crderr.pdb",
        "results/FIKKs_PEXEL/FIKK3/model/model_1.crderr.pdb",
        #expand("results/{prefix}/{sequence_id}/model/model_1.crderr.pdb", sequence_id = sequence_ids),

rule wget_uniprot_fasta:
    output:
        fasta = "results/adhoc/{sequence_id}/{sequence_id}.fasta",
    shell: """
        mkdir -p runs/{wildcards.sequence_id}
        curl https://www.uniprot.org/uniprot/{wildcards.sequence_id}.fasta -o {output.fasta}
    """

rule FIKKs_fasta:
    output:
        fasta = "results/FIKKs/{sequence_id}/{sequence_id}.fasta",
    shell: """
        seqtk subseq results/FIKKs.fasta <(echo "{wildcards.sequence_id}") > {output.fasta}
    """

rule FIKKs_PEXEL_fasta:
    output:
        fasta = "results/FIKKs_PEXEL/{sequence_id}/{sequence_id}.fasta",
    shell: """
        seqtk subseq results/FIKKs_PEXEL.fasta <(echo "{wildcards.sequence_id}") > {output.fasta}
    """

rule run_pyrosetta_ver:
    input:
        fa = "results/{prefix}/{sequence_id}/{sequence_id}.fasta",
    output:
        model_1 = "results/{prefix}/{sequence_id}/model/model_1.crderr.pdb",
    params:
        hhsuite_dir = "/hps/software/users/beltrao/jurgen/RoseTTAFold/software/hhsuite-3.1.0-AVX2-Linux",
    threads: 10
    shell: """
        export PATH="{params.hhsuite_dir}/bin:{params.hhsuite_dir}/scripts:$PATH"
        ./run_pyrosetta_ver.sh {input.fa} results/{wildcards.prefix}/{wildcards.sequence_id}
    """
