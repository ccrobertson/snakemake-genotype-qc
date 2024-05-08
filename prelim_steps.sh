rule make_bed:
    input:
        ped = _raw("Reports.ped"),
        map = _raw("Reports.map"),
        illumina_annot = _resources("illumina_files/infinium-global-diversity-array-8-v1-0_D1.hg19.annotated.txt"),
    output:
        multiext(_results("plink/raw"), ".bed",".bim",".fam"),
    params:
        prefix = _results("plink/raw"),
    shell:
        """
        singularity exec --bind /lab workflow/envs/arushi_general.simg python workflow/scripts/add-dummy-sample-CCR.py --ped {input.ped} --map {input.map} --annot {input.illumina_annot} --prefix {params.prefix}
        plink --ped {params.prefix}.ped --map {params.prefix}.map  --make-bed --out {params.prefix}
        """

rule format:
    input:
        multiext(_results("plink/raw"), ".bed",".bim",".fam"),
        illumina_mapping = _resources("illumina_files/infinium-global-diversity-array-8-v1-0_D1_MappingComment.txt"),
        illumina_rsids = _resources("illumina_files/infinium-global-diversity-array-8-v1-0_D1_b153_rsids.txt"),
        fasta = _resources("hg19/hg19.sorted.fa"),
    output:
        multiext(_results("plink/formatted"), ".bed",".bim",".fam"),
    params:
        outdir = _results("plink"),
        input_prefix = _results("plink/raw"),
        output_prefix = _results("plink/formatted"),
    shell:
        """
        #remove multimapping variants
        awk 'NR>1 && $2!="" {{print $1}}' {input.illumina_mapping} > {params.outdir}/remove.txt
        plink --bfile {params.input_prefix} --exclude {params.outdir}/remove.txt --make-bed --out {params.output_prefix}_tmp

        #update variant ids to dbsnp153
        awk 'BEGIN {{OFS="\t"}} NR>1 && $2!="."' {input.illumina_rsids} | grep -v ',' > {params.outdir}/rename.txt
        plink2 --bfile {params.output_prefix}_tmp --update-name {params.outdir}/rename.txt --make-bed --fa {input.fasta} --out {params.output_prefix}_tmp2

        #remove duplicates
        cut -f 2 {params.output_prefix}_tmp2.bim | sort | uniq -d > {params.outdir}/duplicated.txt
        plink --bfile {params.output_prefix}_tmp2 --exclude {params.outdir}/duplicated.txt --make-bed --out {params.output_prefix}

        #allele freq
        plink --freq --bfile {params.output_prefix} --out {params.output_prefix}
        """

