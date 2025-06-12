rule simulate_phenotype:
  output: "data/pheno/pheno_{h2}_{g}_{n_caus}.pheno", "data/pheno/pheno_{h2}_{g}_{n_caus}.beta"
  container: "docker://uwgac/admix-kit:0.1.2"
  resources:
    mem_mb=5000,
    runtime="30m"
  shell:
    """
    admix simulate-admix-pheno \
      --pfile source_data/CEU-YRI_chrom1 \
      --hsq {wildcards.h2} \
      --cor {wildcards.g} \
      --n-causal {wildcards.n_caus} \
      --out-prefix data/pheno/pheno_{wildcards.h2}_{wildcards.g}_{wildcards.n_caus}
    """

## Add NAs to phenotype file to denote which SNPs to include in PGS
rule fix_beta:
  input: "data/pheno/pheno_{h2}_{g}_{n_caus}.beta"
  output: "data/pheno/pheno_{h2}_{g}_{n_caus}.beta.fixed"
  resources:
    mem_mb=5000,
    runtime="30m"
  shell:
    """
    module load r/4.4.0
    Rscript workflow/scripts/fix_beta.R --input {input} --output {output}
    """

rule train_haudi:
  input:
    pheno="data/pheno/pheno_{h2}_{g}_{n_caus}.pheno",
    beta="data/pheno/pheno_{h2}_{g}_{n_caus}.beta.fixed",
    bk="data/CEU-YRI.bk",
    rds="data/CEU-YRI.rds",
    info="data/CEU-YRI.fbm_info.txt"
  output:
    results=expand("results/haudi/pheno_{{h2}}_{{g}}_{{n_caus}}_{{gamma}}_rep{rep}.txt", rep=range(1, 11)),
    model=expand("model/haudi/pheno_{{h2}}_{{g}}_{{n_caus}}_{{gamma}}_rep{rep}.rds", rep=range(1, 11))
  resources:
    mem_mb=8000,
    runtime="4h"
  shell:
    """
    module load r/4.4.0
    for rep in {{1..10}}; do
      Rscript workflow/scripts/run_haudi.R \
        --gamma {wildcards.gamma} \
        --pheno_file {input.pheno} \
        --beta_file {input.beta} \
        --results_file results/haudi/pheno_{wildcards.h2}_{wildcards.g}_{wildcards.n_caus}_{wildcards.gamma}_rep${{rep}}.txt \
        --model_file model/haudi/pheno_{wildcards.h2}_{wildcards.g}_{wildcards.n_caus}_{wildcards.gamma}_rep${{rep}}.rds \
        --rep ${{rep}}
    done
    """


rule train_gaudi:
  input:
    pheno="data/pheno/pheno_{h2}_{g}_{n_caus}.pheno",
    beta="data/pheno/pheno_{h2}_{g}_{n_caus}.beta.fixed",
    bk="data/CEU-YRI.bk",
    rds="data/CEU-YRI.rds",
    info="data/CEU-YRI.fbm_info.txt"
  output:
    results="results/gaudi/pheno_{h2}_{g}_{n_caus}_{gamma}_rep{rep}.txt",
    model="model/gaudi/pheno_{h2}_{g}_{n_caus}_{gamma}_rep{rep}.rds"
  resources:
    mem_mb=10000,
    runtime="1d"
  shell:
    """
    module load r/4.4.0
    Rscript workflow/scripts/run_gaudi.R \
       --gamma {wildcards.gamma} \
        --pheno_file {input.pheno} \
        --beta_file {input.beta} \
        --results_file {output.results} \
        --model_file {output.model} \
        --rep {wildcards.rep}
    """

rule train_lasso:
  input:
    pheno="data/pheno/pheno_{h2}_{g}_{n_caus}.pheno",
    beta="data/pheno/pheno_{h2}_{g}_{n_caus}.beta.fixed",
    bk="data/CEU-YRI.bk",
    rds="data/CEU-YRI.rds",
    info="data/CEU-YRI.fbm_info.txt"
  output:
    results=expand("results/lasso/pheno_{{h2}}_{{g}}_{{n_caus}}_rep{rep}.txt", rep=range(1, 11)),
    model=expand("model/lasso/pheno_{{h2}}_{{g}}_{{n_caus}}_rep{rep}.rds", rep=range(1, 11))
  resources:
    mem_mb=8000,
    runtime="4h"
  shell:
    """
    module load r/4.4.0
    for rep in {{1..10}}; do
        Rscript workflow/scripts/run_lasso.R \
          --pheno_file {input.pheno} \
          --beta_file {input.beta} \
          --results_file results/lasso/pheno_{wildcards.h2}_{wildcards.g}_{wildcards.n_caus}_rep${{rep}}.txt \
          --model_file model/lasso/pheno_{wildcards.h2}_{wildcards.g}_{wildcards.n_caus}_rep${{rep}}.rds \
          --rep ${{rep}}
    done
    """
