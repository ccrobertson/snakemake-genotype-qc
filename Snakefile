
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


rule all:
    input:
        "results/clean_nosexcheck/clean.bed", 
        

rule sort_chr:
    input:
        config["plink_bfile"] + ".bed"
    output:
        "results/sort_chr/full.bed"
    params:
        input_prefix = config["plink_bfile"],
        output_prefix = "results/sort_chr/full",
    container: "docker://ccrober/r-genotype-qc:latest"
    shell:
        """
        plink1.9 --bfile {params.input_prefix} --allow-extra-chr --make-bed --out {params.output_prefix}
        """

rule missing:
    input:
        "results/sort_chr/full.bed"
    output:
        "results/sort_chr/full.smiss",
        "results/filtered_geno_0.05/tmp.bed",
        "results/filtered_geno_0.05/tmp.smiss",
    params:
        input_prefix = "results/sort_chr/full",
        output_prefix = "results/filtered_geno_0.05/tmp",
    container: "docker://ccrober/r-genotype-qc:latest"
    shell:
        """
        plink2 --bfile {params.input_prefix} --missing --out {params.input_prefix}
        plink2 --bfile {params.input_prefix} --geno 0.05 --mind 0.05 --hwe 1e-6 --make-bed --out {params.output_prefix}
        plink2 --bfile {params.output_prefix} --missing --out {params.output_prefix}
        """

# NOTE: skipping sex check for now because ADI genotypes did not include sex chromosomes
rule sexcheck:
    input: 
        "results/filtered_geno_0.05/tmp.bed"
    output:
        "results/sexcheck/tmp_XY.sexcheck"
    params:
        input_prefix = "results/filtered_geno_0.05/tmp",
        output_prefix = "results/sexcheck/tmp_XY",
    container: "docker://ccrober/r-genotype-qc:latest"
    shell:
        """
        plink1.9 --bfile {params.input_prefix} --allow-extra-chr --chr X,Y --make-bed --out {params.output_prefix}
        plink1.9 --bfile {params.output_prefix} --chr X,Y --check-sex --out {params.output_prefix}
        """

rule pca:
    input: 
        sample_bed = "results/filtered_geno_0.05/tmp.bed",
        ref_bed = "resources/KGref.bed",
    output:
        pcs = "results/pca/kingpcapc.txt",
        popref = "results/pca/kingpca_popref.txt"
    params:
        input_prefix = "results/filtered_geno_0.05/tmp",
        output_prefix = "results/pca/kingpca",
    container: "docker://ccrober/r-genotype-qc:latest"
    shell:
        """
        king -b {input.ref_bed},{input.sample_bed} --pca --projection --prefix {params.output_prefix}
        """

rule ancestry_svm:
    input:
        pcs = "results/pca/kingpcapc.txt",
        popref = "results/pca/kingpca_popref.txt"
    output:
        "results/ancestry_svm/kingsvm_InferredAncestry.txt"
    params:
        output_prefix = "results/ancestry_svm/kingsvm",
    container: "docker://ccrober/r-genotype-qc:latest"
    shell:
        """
        Rscript scripts/Ancestry_Inference.R {input.pcs} {input.popref} {params.output_prefix}
        """

rule relationship_inference:
    input: 
        "results/filtered_geno_0.05/tmp.bed"
    output:
        "results/relationship_inference/kingrel.kin"
    params:
        input_prefix = "results/filtered_geno_0.05/tmp",
        output_prefix = "results/relationship_inference/kingrel",
    container: "docker://ccrober/r-genotype-qc:latest"
    shell:
        """
        #force all donors to have same FID which makes all kinship output save to .kin file
        awk '{{print $1, $2, "0", $2}}' {params.input_prefix}.fam > {params.output_prefix}_update_famids_to_zero.txt
        plink1.9 --bfile {params.input_prefix} --update-ids {params.output_prefix}_update_famids_to_zero.txt --make-bed --out {params.output_prefix}
        king -b {input} --related --prefix {params.output_prefix}
        """

rule summary_nosexcheck:
    input:
        donor_meta_data = config["donor_meta_data"],
        batch_table = config["batch_table"],
        smiss = "results/filtered_geno_0.05/tmp.smiss",
        ancestry_qc = "results/ancestry_svm/kingsvm_InferredAncestry.txt",
        rel_qc = "results/relationship_inference/kingrel.kin",
    output:
        "results/summary_nosexcheck/geno_qc_summary.txt",
        "results/summary_nosexcheck/geno_qc_plots.pdf",
        "results/summary_nosexcheck/duplicates_to_drop.txt",
        "results/summary_nosexcheck/duplicates_to_update.txt",
    params:
        outdir = "results/summary_nosexcheck",
    conda: 
        "rtidy",
    shell: 
        """
        Rscript scripts/qc_summary_nosexcheck.R \
            --donor_meta_data {input.donor_meta_data} \
            --batch_table {input.batch_table} \
            --smiss {input.smiss} \
            --ancestry_qc {input.ancestry_qc} \
            --rel_qc {input.rel_qc} \
            --outdir {params.outdir} 
        """

rule clean_nosexcheck:
    input:
        bed = "results/filtered_geno_0.05/tmp.bed",
        samples_to_drop = "results/summary_nosexcheck/duplicates_to_drop.txt",
        samples_to_update = "results/summary_nosexcheck/duplicates_to_update.txt",
    output:
        "results/clean_nosexcheck/clean.bed", 
    params:
        input_prefix = "results/filtered_geno_0.05/tmp",
        output_prefix = "results/clean_nosexcheck/clean",
    container: "docker://ccrober/r-genotype-qc:latest"
    shell:
        """
        plink1.9 --bfile {params.input_prefix} --remove {input.samples_to_drop} --make-bed --out {params.output_prefix}_TMP1
        plink1.9 --bfile {params.output_prefix}_TMP1 --update-ids {input.samples_to_update} --make-bed --out {params.output_prefix}
        """

### VERSION INCLUDING SEXCHECK
rule summary:
    input:
        donor_meta_data = config["donor_meta_data"],
        batch_table = config["batch_table"],
        smiss = "results/filtered_geno_0.05/tmp.smiss",
        sex_qc = "results/sexcheck/tmp_XY.sexcheck",
        ancestry_qc = "results/ancestry_svm/kingsvm_InferredAncestry.txt",
        rel_qc = "results/relationship_inference/kingrel.kin",
    output:
        "results/summary/geno_qc_summary.txt",
        "results/summary/geno_qc_plots.pdf",
        "results/summary/sex_update.txt",
        "results/summary/duplicates_to_drop.txt",
        "results/summary/duplicates_to_update.txt",
    params:
        outdir = "results/summary",
    conda: "rtidy"
    shell: 
        """
        Rscript scripts/qc_summary.R \
            --donor_meta_data {input.donor_meta_data} \
            --batch_table {input.batch_table} \
            --smiss {input.smiss} \
            --sex_qc {input.sex_qc} \
            --ancestry_qc {input.ancestry_qc} \ 
            --rel_qc {input.rel_qc} \
            --outdir {params.outdir} 
        """

rule clean:
    input:
        bed = "results/filtered_geno_0.05/tmp.bed",
        sex_update = "results/sex_update.txt",
        samples_to_drop = "results/summary/duplicates_to_drop.txt",
        samples_to_update = "results/summary/duplicates_to_update.txt",
    output:
        "results/clean/clean.bed", 
    params:
        input_prefix = "results/filtered_geno_0.05/tmp",
        output_prefix = "results/clean/clean",
    container: "docker://ccrober/r-genotype-qc:latest"
    shell:
        """
        plink1.9 --bfile {params.input_prefix} --update-sex {input.sex_update} --make-bed --out {params.output_prefix}_TMP1
        plink1.9 --bfile {params.output_prefix}_TMP1 --remove {input.samples_to_drop} --make-bed --out {params.output_prefix}_TMP2
        plink1.9 --bfile {params.output_prefix}_TMP2 --update-ids {input.samples_to_update} --make-bed --out {params.output_prefix}
        """

## OLD
# rule summary:
#     input:
#         donor_meta_data = config["donor_meta_data"],
#         batch_table = config["batch_table"],
#         missing_qc = "results/tmp.smiss",
#         sex_qc = "results/tmp_XY.sexcheck",
#         ancestry_qc = "results/kingsvm_InferredAncestry.txt",
#         rel_qc = "results/kingrel.kin",
#     output:
#         "results/geno_qc_summary.txt",
#         "results/geno_qc_plots.pdf",
#         "results/sex_update.txt",
#         "results/duplicates_to_drop.txt",
#         "results/duplicates_to_update.txt",
#     script: "scripts/qc_summary.R"

# rule clean:
#     input:
#         bed = "results/tmp.bed",
#         sex_update = "results/sex_update.txt",
#         samples_to_drop = "results/duplicates_to_drop.txt",
#         samples_to_update = "results/duplicates_to_update.txt",
#     output:
#         "results/clean.bed", 
#     params:
#         input_prefix = "results/tmp",
#         output_prefix = "results/clean",
#     shell:
#         """
#         plink1.9 --bfile {params.input_prefix} --update-sex {input.sex_update} --make-bed --out {params.input_prefix}_TMP1
#         plink1.9 --bfile {params.input_prefix}_TMP1 --remove {input.samples_to_drop} --make-bed --out {params.input_prefix}_TMP2
#         plink1.9 --bfile {params.input_prefix}_TMP2 --update-ids {input.samples_to_update} --make-bed --out {params.output_prefix}
#         """
