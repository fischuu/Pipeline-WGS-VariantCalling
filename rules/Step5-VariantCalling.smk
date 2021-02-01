rule BaseRecalibration:
   """
    Perform the base recalibration (PICARD)
    """
    input:
        bam="%s/BAM/{intid}.sorted.dedup.bam" % (config["project-folder"]),
        bai="%s/BAM/{intid}.sorted.dedup.bam.bai" % (config["project-folder"]),
        ref=config["reference"],
        fai=config["reference-fai"],
        dict=config["reference-dict"],
    output:
        "%s/GATK/recal/{intid}.recal.table" % (config["project-folder"])
    log:
        "%s/logs/GATK/Recalibrate_{intid}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/GATK/Recalibrate_{intid}.benchmark.tsv" % (config["project-folder"])
    params:
        threads=config["params"]["gatk"]["threads"],
        known=config["known-variants"]
    singularity: config["singularity"]["1kbulls"]
    shell:"""
        java -Xmx80G -jar /GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T BaseRecalibrator -nct {params.threads} -R {input.ref} -I {input.bam} -knownSites: {params.known} \
        --bqsrBAQGapOpenPenalty 45 -o {output} &> {log}
    """

rule PrintReads:
   """
    Print reads (PICARD)
    """
    input:
        bam="%s/BAM/{intid}.sorted.dedup.bam" % (config["project-folder"]),
        ref=config["reference"],
        recal="%s/GATK/recal/{intid}.recal.table" % (config["project-folder"])
    output:
        bam="%s/BAM/{intid}.dedub.recal.bam" % (config["project-folder"]),        
        bai="%s/BAM/{intid}.dedub.recal.bam.bai" % (config["project-folder"]),        
    log:
        "%s/logs/GATK/PrintReads_{intid}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/GATK/PrintReads_{intid}.benchmark.tsv" % (config["project-folder"])
    params:
        threads=config["params"]["gatk"]["threads"],
        known=config["known-variants"]
    singularity: config["singularity"]["1kbulls"]
    shell:"""
        java -Xmx80G -jar /GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T PrintReads -nct {params.threads} -R {input.ref} -I {input.bam} -BQSR {input.recal} -o {output.bam}
        
        samtools index {output.bam}
    """
      
rule AnalyzeCovariates:
    """
    Analyze Covariates (PICARD)
    """
    input:
        ref=config["reference"],
        table="%s/GATK/recal/{intid}.recal.table" % (config["project-folder"])
    output:
        table="%s/GATK/recal/{intid}_after_recal.table" % (config["project-folder"]),
        pdf="%s/GATK/recal/{intid}_recal_plots.pdf" % (config["project-folder"])
    log:
        "%s/logs/GATK/AnalyzeCovariates_{intid}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/GATK/AnalyzeCovariates_{intid}.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["1kbulls"]
    shell:"""
        java -Xmx80G -jar /GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T AnalyzeCovariates -R {input.ref} -before {input.table} -after {output.table} -plots {output.pdf} &> {logs}
    """

rule GATK_haplotypeCaller:
    """
    Create the gvcf files (GATK)
    """
    input:
        ref=config["reference"],
        bam="%s/BAM/{intid}.dedub.recal.bam" % (config["project-folder"]),
        bai="%s/BAM/{intid}.dedub.recal.bam.bai" % (config["project-folder"])    
    output:
        "%s/GVCF/{intid}_dedup_recal.g.vcf.gz" % (config["project-folder"])
    log:
        "%s/logs/GATK/HaplotypeCaller_{intid}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/GATK/HaplotypeCaller_{intid}.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["1kbulls"]
    params:
        threads=config["params"]["gatk"]["threads"]
    shell:"""
        java -Xmx80G -jar /GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T HaplotypeCaller -nct {params.threads} -R {input.ref} -I {input.bam} -o {output} -ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000
    """

rule GATK_CallableLoci:
    """
    Get the callable loci (GATK)
    """
    input:
        ref=config["reference"],
        bam="%s/BAM/{intid}.dedub.recal.bam" % (config["project-folder"]),
    output:
        summary="%s/GATK/CallableLoci/{intid}.CallableLoci.summary.txt" % (config["project-folder"]),
        bed="%s/GATK/CallableLoci/{intid}.CallableLoci.bed" % (config["project-folder"])
    log:
        "%s/logs/GATK/CallableLoci_{intid}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/GATK/CallableLoci_{intid}.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["1kbulls"]
    shell:"""
        java -Xmx15g -jar /GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T CallableLoci -R {input.ref} -I {input.bam} -summary {output.summary} -o {output.bed}
    """

rule GATK_DepthOfCoverage:
    """
    Get the depth of coverage (GATK)
    """
    input:
        ref=config["reference"],
        bam="%s/BAM/{intid}.dedub.recal.bam" % (config["project-folder"]),
    output:
        "%s/GATK/DepthOfCoverage/{intid}_dedup_recal.coverage" % (config["project-folder"])
    log:
        "%s/logs/GATK/CallableLoci_{intid}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/GATK/CallableLoci_{intid}.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["1kbulls"]
    shell:"""
        java -Xmx15g -jar  -T CallableLoci -R {input.ref} -I {input.bam} -summary {output.summary} -o {output.bed}
        
        java -Xmx80G -jar /GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T DepthOfCoverage -R {input.ref} -I {input.bam} --omitDepthOutputAtEachBase --logging_level ERROR --summaryCoverageThreshold 10 --summaryCoverageThreshold 20 --summaryCoverageThreshold 30 --summaryCoverageThreshold 40 --summaryCoverageThreshold 50 --summaryCoverageThreshold 80 --summaryCoverageThreshold 90 --summaryCoverageThreshold 100 --summaryCoverageThreshold 150 --minBaseQuality 15 --minMappingQuality 30 --start 1 --stop 1000 --nBins 999 -dt NONE -o {output}
    """
