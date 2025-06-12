rule reference_pgen_to_vcf:
    output:
        vcf="data/CEU-YRI.vcf.gz",
    container: "docker://quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0"
    shell:
        """
        plink2 \
            --pfile source_data/CEU-YRI_chrom1 \
            --keep-allele-order \
            --export vcf 'bgz' \
            --out data/CEU-YRI
        """

rule index_vcf:
    input:
        "data/CEU-YRI.vcf.gz"
    output:
        "data/CEU-YRI.vcf.gz.tbi"
    container: "docker://quay.io/biocontainers/bcftools:1.20--h8b25389_1"
    shell:
        """
        tabix {input}
        """

rule make_switchpoint_vcf:
    input:
        vcf="data/CEU-YRI.vcf.gz",
        index="data/CEU-YRI.vcf.gz.tbi"
    output:
        "data/CEU-YRI.switch.vcf"
    resources:
        mem_mb = 30000,
        runtime = "4h"
    shell:
        """
        module load r/4.4.0
        Rscript workflow/scripts/make_switchpoint_vcf.R \
            --vcf_file {input.vcf} \
            --lanc_file source_data/CEU-YRI_chrom1.lanc \
            --out_file {output}
        """

rule compress_switchpoint_vcf:
    input:
        "data/CEU-YRI.switch.vcf"
    output:
        vcf="data/CEU-YRI.switch.vcf.gz",
        index="data/CEU-YRI.switch.vcf.gz.tbi"
    container: "docker://quay.io/biocontainers/bcftools:1.20--h8b25389_1"
    shell:
        """
        bcftools view {input} -Oz -o {output.vcf}
        tabix {output.vcf}
        """

rule make_lanc_vcf:
    input:
        vcf_switch="data/CEU-YRI.switch.vcf.gz",
        index_switch="data/CEU-YRI.switch.vcf.gz.tbi",
        vcf="data/CEU-YRI.vcf.gz",
        index="data/CEU-YRI.vcf.gz.tbi"
    output:
        vcf="data/CEU-YRI.lanc.vcf.gz",
        index="data/CEU-YRI.lanc.vcf.gz.tbi"
    container: "docker://quay.io/biocontainers/bcftools:1.20--h8b25389_1"
    shell:
        """
        bcftools annotate -c FORMAT -a {input.vcf_switch} {input.vcf} \
            -Oz -o {output.vcf}
        tabix {output.vcf}
        """

rule make_fbm:
    input:
        "data/CEU-YRI.lanc.vcf.gz",
        "data/CEU-YRI.lanc.vcf.gz.tbi"
    output:
        "data/CEU-YRI.bk",
        "data/CEU-YRI.rds",
        "data/CEU-YRI.fbm_info.txt"
    resources:
        mem_mb=10000,
        runtime="1h"
    shell:
        """
        module load r/4.4.0
        Rscript workflow/scripts/make_FBM.R
        """
