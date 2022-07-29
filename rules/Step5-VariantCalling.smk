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
        threads=config["params"]["gatk"]["threads"],
        known=config["known-variants"]
    singularity: config["singularity"]["wgs"]
    shell:"""
        java -Xmx80G -jar /GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T BaseRecalibrator -nct {params.threads} -R {input.ref} -I {input.bam} -knownSites: {params.known} \
        --bqsrBAQGapOpenPenalty 45 -o {output} &> {log}
    """

rule PrintReads:
   """
    Print reads (PICARD)
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
        "%s/logs/GATK/PrintReads_{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/GATK/PrintReads_{samples}.benchmark.tsv" % (config["project-folder"])
    params:
        threads=config["params"]["gatk"]["threads"],
        known=config["known-variants"]
    singularity: config["singularity"]["wgs"]
    shell:"""
        gatk --java-options '-Xmx80G' PrintReads -nct {params.threads} -R {input.ref} -I {input.bam} -BQSR {input.recal} -o {output.bam} &> {log}
        
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
        threads=config["params"]["gatk"]["threads"],
        known=config["known-variants"]
    singularity: config["singularity"]["wgs"]
    shell:"""
        gatk --java-options '-Xmx80G' BaseRecalibrator -nct {params.threads} \
            -R {input.ref} \
            -I {input.bam} \
            -knownSites: {params.known} \
            --bqsrBAQGapOpenPenalty 45 \
            -o {output} &> {log}
    """

      
rule AnalyzeCovariates:
    """
    Analyze Covariates (PICARD)
    """
    input:
        ref=config["reference"],
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
             -R {input.ref} \
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
        gatk --java-options '-Xmx80G' HaplotypeCaller -nct {params.threads} \
             -R {input.ref} \
             -I {input.bam} \
             -o {output.vcf} \
             -ERC GVCF \
             -variant_index_type LINEAR \
             -variant_index_parameter 128000
        
        md5sum {output.vcf} > {output.md5}
    """

rule GATK_CallableLoci:
    """
    Get the callable loci (GATK)
    """
    input:
        ref=config["reference"],
        bam="%s/BAM/{samples}.dedup.recal.bam" % (config["project-folder"]),
    output:
        summary="%s/GATK/CallableLoci/{samples}.CallableLoci.summary.txt" % (config["project-folder"]),
        bed="%s/GATK/CallableLoci/{samples}.CallableLoci.bed" % (config["project-folder"]),
        summarymd5="%s/GATK/CallableLoci/{samples}.CallableLoci.summary.txt.md5" % (config["project-folder"]),
        bedmd5="%s/GATK/CallableLoci/{samples}.CallableLoci.bed.md5" % (config["project-folder"])
    log:
        "%s/logs/GATK/CallableLoci_{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/GATK/CallableLoci_{samples}.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["wgs"]
    shell:"""
        gatk --java-options '-Xmx15g' CallableLoci \
             -R {input.ref} \
             -I {input.bam} \
             -summary {output.summary} \
             -o {output.bed}
        
        md5sum {output.summary} > {output.summarymd5}
        md5sum {output.bed} > {output.bedmd5}
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
        "%s/logs/GATK/CallableLoci_{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/GATK/CallableLoci_{samples}.benchmark.tsv" % (config["project-folder"])
    params: out="%s/GATK/DepthOfCoverage/{samples}_dedup_recal.coverage" % (config["project-folder"])
    singularity: config["singularity"]["wgs"]
    shell:"""
        gatk --java-options '-Xmx80G' DepthOfCoverage \
             -R {input.ref} \
             -I {input.bam} \
             --omitDepthOutputAtEachBase \
             --logging_level ERROR \
             --summaryCoverageThreshold 10 \
             --summaryCoverageThreshold 20 \
             --summaryCoverageThreshold 30 \
             --summaryCoverageThreshold 40 \
             --summaryCoverageThreshold 50 \
             --summaryCoverageThreshold 80 \
             --summaryCoverageThreshold 90 \
             --summaryCoverageThreshold 100 \
             --summaryCoverageThreshold 150 \
             --minBaseQuality 15 \
             --minMappingQuality 30 \
             --start 1 \
             --stop 1000 \
             --nBins 999 \
             -dt NONE \
             -o {params.out} &> {log}
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
             --out {output.gvcf} &> {log}
        
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
             --out {output.vcf} &> {log}
        
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
             --resource:knownvariants,known=true,training=true,truth=false,prior=15 {input.variants} \
             -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
             --out {output.reca} \
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
    log:
        "%s/logs/GATK/VariantRecalibrator.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/GATK/VariantRecalibrator.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["wgs"]
    shell:"""
        gatk --java-options '-Xmx80G' ApplyVQSR \
             -R {input.ref} \
             -V {input.vcf} \
             --out {output.vcf} \
             --truth-sensitivity-filter-level 99.0 \
             --tranches-file {input.tranches} \
             --recal-file {input.reca} \
             -mode SNP &> {log}
    """