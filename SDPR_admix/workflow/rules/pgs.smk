rule pgs:
  input:
    vcf_train = "data/geno_subset/chr1_train_{HSQ}_{GEN_COR}_{N_CAUS}_rep{REP}.vcf.gz",
    vcf_test = "data/geno_subset/chr1_test_{HSQ}_{GEN_COR}_{N_CAUS}_rep{REP}.vcf.gz",
    msp_train = "data/geno_subset/chr1_train_{HSQ}_{GEN_COR}_{N_CAUS}_rep{REP}.msp.tsv",
    msp_test = "data/geno_subset/chr1_test_{HSQ}_{GEN_COR}_{N_CAUS}_rep{REP}.msp.tsv",
    pheno = "data/pheno/train/train_{HSQ}_{GEN_COR}_{N_CAUS}_rep{REP}.pheno"
  output:
    model = "model/result_{HSQ}_{GEN_COR}_{N_CAUS}_rep{REP}.txt",
    score = "score/score_{HSQ}_{GEN_COR}_{N_CAUS}_rep{REP}.txt"
  resources:
    mem_mb = 10000,
    runtime = "1h"
  params:
    rho=lambda wildcards: "0.99" if float(wildcards.GEN_COR) == 1 else wildcards.GEN_COR
  container: "docker://frankpo/run_sdpr_admix:0.0.1"
  shell:
    """
    SDPR_admix \
      -vcf {input.vcf_train} \
      -pheno {input.pheno} \
      -msp {input.msp_train} \
      -rho {params.rho} \
      -out {output.model}

    score \
      -vcf {input.vcf_test} \
      -msp {input.msp_test} \
      -score {output.model} \
      -out {output.score}
    """
