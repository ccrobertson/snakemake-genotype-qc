
################################################
#
# Provide config file with the following format:
#
# donor_phenotypes: <path_to_pheno_file> (must include the variables "Donor" and "Sex")>
# plink_prefix: <path_to_plink_files>
#
# plink_bfile: <path_to_plink_files>
# donor_meta_data: <path_to_pheno_file> (must include the variables "Donor" and "Sex")>
# batch_table: <path_to_genotyping_batch_files> (must include IID, Donor, Batch, chip_version)
#
################################################
configfile: "config.yaml"
container: "docker://ccrober/r-genotype-qc:latest"

rule all:
    input:
        #"results/full.bed", 
        #"results/tmp_XY.sexcheck",
        #"results/kingpcapc.txt",
        #"results/kingsvm_InferredAncestry.txt",
        #"results/kingrel.kin",
        #"results/geno_qc_summary.txt",
        #"results/geno_qc_plots.pdf",
        "results/clean.bed"
        

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
        pcs = "results/kingpcapc.txt",
        popref = "results/kingpca_popref.txt"
    shell:
        """
        king -b resources/KGref.bed,{input} --pca --projection --prefix "results/kingpca"
        """

rule ancestry_svm:
    input:
        pcs = "results/kingpcapc.txt",
        popref = "results/kingpca_popref.txt"
    output:
        "results/kingsvm_InferredAncestry.txt"
    shell:
        """
        Rscript scripts/Ancestry_Inference.R {input.pcs} {input.popref} results/kingsvm
        """

rule related:
    input: 
        "results/tmp.bed"
    output:
        "results/kingrel.kin"
    shell:
        """
        king -b {input} --related --prefix "results/kingrel"
        """

rule summary:
    input:
        donor_meta_data = config["donor_meta_data"],
        batch_table = config["batch_table"],
        missing_qc = "results/tmp.smiss",
        sex_qc = "results/tmp_XY.sexcheck",
        ancestry_qc = "results/kingsvm_InferredAncestry.txt",
        rel_qc = "results/kingrel.kin",
    output:
        "results/geno_qc_summary.txt",
        "results/geno_qc_plots.pdf",
        "results/sex_update.txt",
        "results/duplicates_to_drop.txt",
        "results/duplicates_to_update.txt",
    script: "scripts/qc_summary.R"

rule clean:
    input:
        bed = "results/tmp.bed",
        sex_update = "results/sex_update.txt",
        samples_to_drop = "results/duplicates_to_drop.txt",
        samples_to_update = "results/duplicates_to_update.txt",
    output:
        "results/clean.bed", 
    params:
        input_prefix = "results/tmp",
        output_prefix = "results/clean",
    shell:
        """
        plink1.9 --bfile {params.input_prefix} --update-sex {input.sex_update} --make-bed --out {params.input_prefix}_TMP1
        plink1.9 --bfile {params.input_prefix}_TMP1 --remove {input.samples_to_drop} --make-bed --out {params.input_prefix}_TMP2
        plink1.9 --bfile {params.input_prefix}_TMP2 --update-ids {input.samples_to_update} --make-bed --out {params.output_prefix}
        """
