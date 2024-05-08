
################################################
#
# Provide config file with the following format:
#
# donor_phenotypes: <path_to_pheno_file> (must include the variables "Donor" and "Sex")>
# plink_prefix: <path_to_plink_files>
#
################################################
configfile: "config.yaml"

rule all:
    input: 
        "results/tmp_XY.sexcheck",
        "results/tmp_kingpc.txt",
        "results/tmp_king.kin",
        "results/tmp_king_InferredAncestry.txt",

rule sort_chr:
    input:
        config["plink_bfile"] + ".bed"
    output:
        "results/full.bed"
    params:
        prefix = config["plink_bfile"], 
    shell:
        """
        plink1.9 --bfile {params.prefix} --allow-extra-chr --make-bed --out results/full
        """

rule missing:
    input:
        "results/full.bed"
    output:
        "results/full.smiss",
        "results/tmp.bed",
        "results/tmp.smiss",
    shell:
        """
        plink2 --bfile results/full --missing --out results/full
        plink2 --bfile results/full --geno 0.05 --mind 0.05 --hwe 1e-6 --make-bed --out results/tmp
        plink2 --bfile results/tmp --missing --out results/tmp
        """

rule sexcheck:
    input: 
        "results/tmp.bed"
    output:
        "results/tmp_XY.sexcheck"
    shell:
        """
        plink1.9 --bfile results/tmp --allow-extra-chr --chr X,Y --make-bed --out results/tmp_XY
        plink1.9 --bfile results/tmp_XY --chr X,Y --check-sex --out results/tmp_XY
        """

rule pca:
    input: 
        "results/tmp.bed"
    output:
        "results/tmp_kingpc.txt"
    shell:
        """
        king -b resources/KGref.bed,{input} --pca --projection --prefix "results/tmp_king"
        """

rule ancestry_svm:
    input:
        pcs = "results/tmp_kingpc.txt",
        popref = "results/tmp_king_popref.txt"
    output:
        "results/tmp_king_InferredAncestry.txt"
    conda:
        "renv"
    shell:
        """
        Rscript scripts/Ancestry_Inference.R {input.pcs} {input.popref} results/tmp_king
        """

rule related:
    input: 
        "results/tmp.bed"
    output:
        "results/tmp_king.kin"
    shell:
        """
        king -b {input} --related --prefix "results/tmp_king"
        """

## NEXT MAKE REPORT TO SUMMARIZE QC OUTPUT
# 1. SEX CHECK OVERVIEW
# 2. ANCESTRY INFERRED VS REPORTED
# 3. ANY RELATEDNESS?
# rule report:
#     input:
#         pheno = config["donor_phenotypes"],
#         pcs = "results/kingpc.txt"
#     output:
#     conda:
#         "renv"
#     shell:
#         """
#         """
