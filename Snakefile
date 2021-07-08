
rule all:
    input:
        "runs/O00244/model/model_1.crderr.pdb",
        "runs/P42081/model/model_1.crderr.pdb",
        "runs/Q13148/model/model_1.crderr.pdb",

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
