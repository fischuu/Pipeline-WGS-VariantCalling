rule BaseRecalibration:
   """
    Perform the base recalibration (PICARD)
    """
    input:
        bam="%s/BAM/{samples}.sorted.dedup.bam" % (config["project-folder"]),
        bai="%s/BAM/{samples}.sorted.dedup.bam.bai" % (config["project-folder"]),
        ref=config["reference"],
        fai=config["reference-fai"],
        dict=config["reference-dict"],
    output:
        "%s/GATK/recal/{samples}.recal.table" % (config["project-folder"])
    log:
        "%s/logs/GATK/Recalibrate_{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/GATK/Recalibrate_{samples}.benchmark.tsv" % (config["project-folder"])
    params:
        bqsrBAQGapOpenPenalty=config["params"]["gatk"]["bqsrBAQGapOpenPenalty"],
        known=config["known-variants"]
    singularity: config["singularity"]["wgs"]
    shell:"""
        gatk --java-options '-Xmx80G' BaseRecalibrator \
             -R {input.ref} \
             -I {input.bam} \
             --known-sites {params.known} \
             --bqsr-baq-gap-open-penalty {params.bqsrBAQGapOpenPenalty} \
             --output {output} &> {log}
    """

rule ApplyBQRS:
   """
    Apply the model to adjust the base quality scores (PICARD)
    """
    input:
        bam="%s/BAM/{samples}.sorted.dedup.bam" % (config["project-folder"]),
        bai="%s/BAM/{samples}.sorted.dedup.bam.bai" % (config["project-folder"]),
        ref=config["reference"],
        recal="%s/GATK/recal/{samples}.recal.table" % (config["project-folder"])
    output:
        bam="%s/BAM/{samples}.dedup.recal.bam" % (config["project-folder"]),        
        bai="%s/BAM/{samples}.dedup.recal.bam.bai" % (config["project-folder"]),
        md5="%s/BAM/{samples}.dedup.recal.bam.md5" % (config["project-folder"])
    log:
        "%s/logs/GATK/ApplyBQRS_{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/GATK/ApplyBQRS_{samples}.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["wgs"]
    shell:"""
        gatk --java-options '-Xmx80G' ApplyBQSR \
             -I {input.bam} \
             -R {input.ref} \
             --bqsr-recal-file {input.recal} \
             -O {output.bam} &> {log}
    
        samtools index {output.bam}
        md5sum {output.bam} > {output.md5}
    """


rule BaseRecalibration_afterRecal:
   """
    Perform the base recalibration (PICARD)
    """
    input:
        bam="%s/BAM/{samples}.dedup.recal.bam" % (config["project-folder"]),
        bai="%s/BAM/{samples}.dedup.recal.bam.bai" % (config["project-folder"]),
        ref=config["reference"],
        fai=config["reference-fai"],
        dict=config["reference-dict"],
    output:
        "%s/GATK/recal/{samples}_after_recal.table" % (config["project-folder"])
    log:
        "%s/logs/GATK/Recalibrate_after_{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/GATK/Recalibrate_after_{samples}.benchmark.tsv" % (config["project-folder"])
    params:
        bqsrBAQGapOpenPenalty=config["params"]["gatk"]["bqsrBAQGapOpenPenalty"],
        known=config["known-variants"]
    singularity: config["singularity"]["wgs"]
    shell:"""
        gatk --java-options '-Xmx80G' BaseRecalibrator \
            -R {input.ref} \
            -I {input.bam} \
            --known-sites {params.known} \
            --bqsr-baq-gap-open-penalty {params.bqsrBAQGapOpenPenalty} \
            --output {output} &> {log}
    """

      
rule AnalyzeCovariates:
    """
    Analyze Covariates (PICARD)
    """
    input:
        tableBefore="%s/GATK/recal/{samples}.recal.table" % (config["project-folder"]),
        tableAfter="%s/GATK/recal/{samples}_after_recal.table" % (config["project-folder"])
    output:
        pdf="%s/GATK/recal/{samples}_recal_plots.pdf" % (config["project-folder"])
    log:
        "%s/logs/GATK/AnalyzeCovariates_{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/GATK/AnalyzeCovariates_{samples}.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["wgs"]
    shell:"""
        gatk --java-options '-Xmx80G' AnalyzeCovariates \
             -before {input.tableBefore} \
             -after {input.tableAfter} \
             -plots {output.pdf} &> {log}
    """

rule GATK_haplotypeCaller:
    """
    Create the gvcf files (GATK)
    """
    input:
        ref=config["reference"],
        bam="%s/BAM/{samples}.dedup.recal.bam" % (config["project-folder"]),
        bai="%s/BAM/{samples}.dedup.recal.bam.bai" % (config["project-folder"])    
    output:
        vcf="%s/GATK/GVCF/{samples}_dedup_recal.g.vcf.gz" % (config["project-folder"]),
        md5="%s/GATK/GVCF/{samples}_dedup_recal.g.vcf.gz.md5" % (config["project-folder"])
    log:
        "%s/logs/GATK/HaplotypeCaller_{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/GATK/HaplotypeCaller_{samples}.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["wgs"]
    params:
        threads=config["params"]["gatk"]["threads"]
    shell:"""
        gatk --java-options '-Xmx80G' HaplotypeCaller \
             -R {input.ref} \
             -I {input.bam} \
             --output {output.vcf} \
             -ERC GVCF &> {log}
        
        md5sum {output.vcf} > {output.md5}
    """

rule GATK_DepthOfCoverage:
    """
    Get the depth of coverage (GATK)
    """
    input:
        ref=config["reference"],
        bam="%s/BAM/{samples}.dedup.recal.bam" % (config["project-folder"]),
    output:
        "%s/GATK/DepthOfCoverage/{samples}_dedup_recal.coverage.sample_summary" % (config["project-folder"])
    log:
        "%s/logs/GATK/DepthOfCoverage_{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/GATK/DepthOfCoverage_{samples}.benchmark.tsv" % (config["project-folder"])
    params: 
        out="%s/GATK/DepthOfCoverage/{samples}_dedup_recal.coverage" % (config["project-folder"]),
        chr=config["params"]["gatk"]["DepthOfCoverage"]["interval"]
    singularity: config["singularity"]["wgs"]
    shell:"""
        gatk --java-options '-Xmx80G' DepthOfCoverage \
             -R {input.ref} \
             -I {input.bam} \
             -L {params.chr} \
             --omit-depth-output-at-each-base \
             --verbosity ERROR \
             --summary-coverage-threshold 10 \
             --summary-coverage-threshold 20 \
             --summary-coverage-threshold 30 \
             --summary-coverage-threshold 40 \
             --summary-coverage-threshold 50 \
             --summary-coverage-threshold 80 \
             --summary-coverage-threshold 90 \
             --summary-coverage-threshold 100 \
             --summary-coverage-threshold 150 \
             --min-base-quality 15 \
             --start 1 \
             --stop 1000 \
             --nBins 999 \
             --output {params.out} &> {log}
    """

rule GATK_combineGVCFs:
    """
    Get the depth of coverage (GATK)
    """
    input:
        ref=config["reference"],
        gvcfs=expand("%s/GATK/GVCF/{samples}_dedup_recal.g.vcf.gz" % (config["project-folder"]), samples=samples)
    output:
        gvcf="%s/GATK/Cohort.g.vcf.gz" % (config["project-folder"]),
        md5="%s/GATK/Cohort.g.vcf.gz.md5" % (config["project-folder"])
    log:
        "%s/logs/GATK/combineGVCFs.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/GATK/combineGVCFs.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["wgs"]
    shell:"""
        gatk --java-options '-Xmx80G' CombineGVCFs \
             -R {input.ref} \
             --variant $(echo {input.gvcfs} | sed 's/ / --variant /g') \
             --output {output.gvcf} &> {log}
        
        md5sum {output.gvcf} > {output.md5}
    """
    
    
rule GATK_GenotypeGVCFs:
    """
    Perform the genotyping (GATK)
    """
    input:
        ref=config["reference"],
        gvcf="%s/GATK/Cohort.g.vcf.gz" % (config["project-folder"])
    output:
        vcf="%s/GATK/output.vcf.gz" % (config["project-folder"]),
        md5="%s/GATK/output.vcf.gz.md5" % (config["project-folder"])
    log:
        "%s/logs/GATK/genotypeGVCFs.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/GATK/genotypeGVCFs.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["wgs"]
    shell:"""
        gatk --java-options '-Xmx80G' GenotypeGVCFs \
             -R {input.ref} \
             -V {input.gvcf} \
             --output {output.vcf} &> {log}
        
        md5sum {output.vcf} > {output.md5}
    """
    
rule GATK_VariantRecalibrator:
    """
    Build a recalibration model to score variant quality for filtering purposes (GATK)
    """
    input:
        vcf="%s/GATK/output.vcf.gz" % (config["project-folder"]),
        variants=config["known-variants"]
    output:
        reca="%s/GATK/cohort_snps.reca" % (config["project-folder"]),
        tranches="%s/GATK/cohort_snps.tranches" % (config["project-folder"])
    log:
        "%s/logs/GATK/VariantRecalibrator.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/GATK/VariantRecalibrator.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["wgs"]
    shell:"""
        gatk --java-options '-Xmx80G' VariantRecalibrator \
             -V {input.vcf} \
             -mode SNP \
             --max-gaussians 8 \
             --resource:knownvariants,known=false,training=true,truth=true,prior=15 {input.variants} \
             -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
             --output {output.reca} \
             --tranches-file {output.tranches} &> {log}
    """

rule GATK_ApplyVQSR:
    """
    Apply a score cutoff to filter variants based on a recalibration table (GATK)
    """
    input:
        ref=config["reference"],
        vcf="%s/GATK/output.vcf.gz" % (config["project-folder"]),
        reca="%s/GATK/cohort_snps.reca" % (config["project-folder"]),
        tranches="%s/GATK/cohort_snps.tranches" % (config["project-folder"])
    output:
        vcf="%s/GATK/output.vqsr.vcf" % (config["project-folder"]),
        res="%s/RESULTS/final_variants.vcf" % (config["project-folder"])
    log:
        "%s/logs/GATK/VariantRecalibrator.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/GATK/VariantRecalibrator.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["wgs"]
    shell:"""
        gatk --java-options '-Xmx80G' ApplyVQSR \
             -R {input.ref} \
             -V {input.vcf} \
             --output {output.vcf} \
             --truth-sensitivity-filter-level 99.0 \
             --tranches-file {input.tranches} \
             --recal-file {input.reca} \
             -mode SNP &> {log}
             
        cp {output.vcf} {output.res} 
    """