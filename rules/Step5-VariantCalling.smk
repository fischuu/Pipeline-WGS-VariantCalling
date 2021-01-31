rule BaseRecalibration:
   """
    Perform the base recalibration (PICARD)
    """
    input:
        bam="%s/BAM/{intid}.sorted.dedup.bam" % (config["project-folder"]),
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
        "%s/BAM/{intid}.dedub.recal.bam" % (config["project-folder"])        
    log:
        "%s/logs/GATK/PrintReads_{intid}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/GATK/PrintReads_{intid}.benchmark.tsv" % (config["project-folder"])
    params:
        threads=config["params"]["gatk"]["threads"],
        known=config["known-variants"]
    singularity: config["singularity"]["1kbulls"]
    shell:"""
        java -Xmx80G -jar GATK.jar -T PrintReads -nct {params.threads} -R {input.ref} -I {input.bam} -BQSR {input.recal} -o {output}
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
        java -Xmx80G -jar $GATK.jar -T AnalyzeCovariates -R {input.ref} -before {input.table} -after {output.table} -plots {output.pdf}
    """
