# snakemake --cores 10 --dry-run
# snakemake --cores 10 --until wget_uniprot_fasta --dry-run
# ./smk_codon --dry-run

import pandas as pd

def get_sequence_ids(fp):
    return [sequence_id for sequence_id in shell("""grep "^>" %s | awk 'sub(/^>/, "")'""" % (fp,), iterable=True)]

rule all:
    input:
        expand("results/FIKKs/{sequence_id}/model/model_1.crderr.pdb", sequence_id = get_sequence_ids("results/FIKKs.fasta")),
        expand("results/FIKKs_PEXEL/{sequence_id}/model/model_1.crderr.pdb", sequence_id = get_sequence_ids("results/FIKKs_PEXEL.fasta")),

rule wget_uniprot_fasta:
    output:
        fasta = "results/adhoc/{sequence_id}/{sequence_id}.fasta",
    group: "run_pyrosetta_ver"
    shell: """
        mkdir -p runs/{wildcards.sequence_id}
        curl https://www.uniprot.org/uniprot/{wildcards.sequence_id}.fasta -o {output.fasta}
    """

rule FIKKs_fasta:
    output:
        fasta = "results/FIKKs/{sequence_id}/{sequence_id}.fasta",
    group: "run_pyrosetta_ver"
    shell: """
        seqtk subseq results/FIKKs.fasta <(echo "{wildcards.sequence_id}") > {output.fasta}
    """

rule FIKKs_PEXEL_fasta:
    output:
        fasta = "results/FIKKs_PEXEL/{sequence_id}/{sequence_id}.fasta",
    group: "run_pyrosetta_ver"
    shell: """
        seqtk subseq results/FIKKs_PEXEL.fasta <(echo "{wildcards.sequence_id}") > {output.fasta}
    """

rule run_pyrosetta_ver:
    input:
        fa = "results/{prefix}/{sequence_id}/{sequence_id}.fasta",
    output:
        hhblits = directory("results/{prefix}/{sequence_id}/hhblits"),
        model = directory("results/{prefix}/{sequence_id}/model"),
        model1 = "results/{prefix}/{sequence_id}/model/model_1.crderr.pdb",
        model2 = "results/{prefix}/{sequence_id}/model/model_2.crderr.pdb",
        model3 = "results/{prefix}/{sequence_id}/model/model_3.crderr.pdb",
        model4 = "results/{prefix}/{sequence_id}/model/model_4.crderr.pdb",
        model5 = "results/{prefix}/{sequence_id}/model/model_5.crderr.pdb",
        pdb_3track = directory("results/{prefix}/{sequence_id}/pdb-3track"),
        _fold_list = "results/{prefix}/{sequence_id}/parallel.fold.list",
        _3track_npz = "results/{prefix}/{sequence_id}/t000_.3track.npz",
        _atab = "results/{prefix}/{sequence_id}/t000_.atab",
        _hhr = "results/{prefix}/{sequence_id}/t000_.hhr",
        _msa0_a3m = "results/{prefix}/{sequence_id}/t000_.msa0.a3m",
        _ss2_a3m = "results/{prefix}/{sequence_id}/t000_.msa0.ss2.a3m",
        _ss2 = "results/{prefix}/{sequence_id}/t000_.ss2",
    params:
        hhsuite_dir = "/hps/software/users/beltrao/jurgen/RoseTTAFold/software/hhsuite-3.1.0-AVX2-Linux",
    group: "run_pyrosetta_ver"
    threads: 10
    shell: """
        export PATH="{params.hhsuite_dir}/bin:{params.hhsuite_dir}/scripts:$PATH"
        ./run_pyrosetta_ver.sh {input.fa} results/{wildcards.prefix}/{wildcards.sequence_id}
    """
