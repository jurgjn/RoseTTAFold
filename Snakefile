# snakemake --cores 10 --dry-run
# snakemake --cores 10 --until wget_uniprot_fasta --dry-run

import pandas as pd

sequence_ids = pd.read_csv("runs/unstructured_dms_tabfix.tsv", sep='\t')['Uniprot ID']

rule all:
    input:
        expand("runs/{sequence_id}/model/model_1.crderr.pdb", sequence_id = sequence_ids),

rule wget_uniprot_fasta:
    output:
        fasta = "runs/{sequence_id}/{sequence_id}.fasta",
    shell: """
        mkdir -p runs/{wildcards.sequence_id}
        curl https://www.uniprot.org/uniprot/{wildcards.sequence_id}.fasta -o {output.fasta}
    """

rule run_pyrosetta_ver:
    input:
        fa = "runs/{sequence_id}/{sequence_id}.fasta",
    output:
        model_1 = "runs/{sequence_id}/model/model_1.crderr.pdb",
    params:
        hhsuite_dir = "/hps/nobackup/beltrao/jurgen/nobackup/hhsuite-3.1.0-AVX2-Linux",
    threads: 10
    shell: """
        export PATH="{params.hhsuite_dir}/bin:{params.hhsuite_dir}/scripts:$PATH"
        ./run_pyrosetta_ver.sh {input.fa} runs/{wildcards.sequence_id}
    """
