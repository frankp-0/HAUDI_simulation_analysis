rule make_geno_input:
  output:
    vcf = "data/chr1.vcf.gz",
    msp = "data/chr1.msp.tsv"
  resources:
    mem_mb = 30000,
    runtime = "1h"
  params:
    map_file = "/work/users/f/r/frocko/HAUDI/simulation/data/1kg-ref-hg38/metadata/genetic_map/chr1.map",
    plink_prefix = "/proj/yunligrp/users/frank/projects/HAUDI/HAUDI_simulation/source_data/CEU-YRI_chrom1"
  shell:
    """
    plink2 --pfile {params.plink_prefix} --recode vcf bgz --out data/chr1
    Rscript workflow/scripts/make_msp.R --plink_prefix {params.plink_prefix} --map_file {params.map_file} --msp_out {output.msp}
    """

rule list_samples:
  output:
    expand("data/samples/samples_{SPLIT}_rep{REP}.txt", SPLIT = ["train", "test"], REP = REP)
  shell:
     """
     Rscript workflow/scripts/list_samples.R
     """

rule make_pheno:
  output:
    expand(
      "data/pheno/{SPLIT}/{SPLIT}_{HSQ}_{GEN_COR}_{N_CAUS}_rep{REP}.pheno",
      HSQ=HSQ, GEN_COR=GEN_COR, N_CAUS=N_CAUS, REP=REP, SPLIT = ["train", "test"])
  resources:
    mem_mb = 4000,
    runtime = "10m"
  shell:
    """
    Rscript workflow/scripts/make_pheno.R
    """

rule make_subset_geno:
  input:
    vcf = "data/chr1.vcf.gz",
    msp = "data/chr1.msp.tsv",
    samples = expand("data/samples/samples_{SPLIT}_rep{{REP}}.txt", SPLIT = ["train", "test"])
  output:
    vcf = expand("data/geno_subset/chr1_{SPLIT}_{{HSQ}}_{{GEN_COR}}_{{N_CAUS}}_rep{{REP}}.vcf.gz", SPLIT = ["train", "test"]),
    plink = temp(expand("data/geno_subset/chr1_{{HSQ}}_{{GEN_COR}}_{{N_CAUS}}_rep{{REP}}.{PLINK}", PLINK = ["pgen", "psam", "pvar", "log"], SPLIT = ["train", "test"])),
    snps = "data/geno_subset/chr1_{HSQ}_{GEN_COR}_{N_CAUS}_rep{REP}.snps",
    msp = expand("data/geno_subset/chr1_{SPLIT}_{{HSQ}}_{{GEN_COR}}_{{N_CAUS}}_rep{{REP}}.msp.tsv", SPLIT=["train", "test"])
  resources:
    mem_mb = 10000,
    runtime = "20m"
  shell:
    """
    Rscript workflow/scripts/get_pheno_snps.R \
      --beta ../data/pheno/pheno_{wildcards.HSQ}_{wildcards.GEN_COR}_{wildcards.N_CAUS}.beta.fixed \
      --rep {wildcards.REP} \
      --out {output.snps}
    plink2 \
      --vcf {input.vcf} --extract {output.snps} --make-pgen \
      --out data/geno_subset/chr1_{wildcards.HSQ}_{wildcards.GEN_COR}_{wildcards.N_CAUS}_rep{wildcards.REP}
    plink2 --pfile data/geno_subset/chr1_{wildcards.HSQ}_{wildcards.GEN_COR}_{wildcards.N_CAUS}_rep{wildcards.REP} \
      --recode vcf bgz --keep data/samples/samples_train_rep{wildcards.REP}.txt \
      --out data/geno_subset/chr1_train_{wildcards.HSQ}_{wildcards.GEN_COR}_{wildcards.N_CAUS}_rep{wildcards.REP}
    plink2 --pfile data/geno_subset/chr1_{wildcards.HSQ}_{wildcards.GEN_COR}_{wildcards.N_CAUS}_rep{wildcards.REP} \
      --recode vcf bgz --keep data/samples/samples_test_rep{wildcards.REP}.txt \
      --out data/geno_subset/chr1_test_{wildcards.HSQ}_{wildcards.GEN_COR}_{wildcards.N_CAUS}_rep{wildcards.REP}
    Rscript workflow/scripts/reduce_msp.R --msp_in {input.msp} --snps {output.snps} \
      --samples data/samples/samples_train_rep{wildcards.REP}.txt \
      --msp_out data/geno_subset/chr1_train_{wildcards.HSQ}_{wildcards.GEN_COR}_{wildcards.N_CAUS}_rep{wildcards.REP}.msp.tsv
    Rscript workflow/scripts/reduce_msp.R --msp_in {input.msp} --snps {output.snps} \
      --samples data/samples/samples_test_rep{wildcards.REP}.txt \
      --msp_out data/geno_subset/chr1_test_{wildcards.HSQ}_{wildcards.GEN_COR}_{wildcards.N_CAUS}_rep{wildcards.REP}.msp.tsv
    """
