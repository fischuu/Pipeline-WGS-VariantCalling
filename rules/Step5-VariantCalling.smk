rule BaseRecalibration:
   """
    Perform the base recalibration (PICARD)
    """
    input:
        bam="%s/BAM/"${OutputFile}-pe.dedup.bam" % (config["project-folder"]),
        ref=config["reference"]
    output:
        "%s/GATK/recal/${OutputFile}.recal.table" % (config["project-folder"])
    log:
        "%s/logs/GATK/Recalibrate_{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/GATK/Recalibrate_{samples}.benchmark.tsv" % (config["project-folder"])
    params:
        threads=config["params"]["gatk"]["threads"]
        known=config["known-variants"]
    singularity: config["singularity"]["1kbulls"]
    shell:"""
        java -Xmx80G -jar /path/to/GATK.jar –T BaseRecalibrator –nct {params.threads} -R {input.ref} -I {input.bam} –knownSites:vcf {params.known} \
        -–bqsrBAQGapOpenPenalty 45 -o {output}
    """

rule PrintReads:
   """
    Print reads (PICARD)
    """
    input:
        bam="%s/BAM/"${OutputFile}-pe.dedub.bam" % (config["project-folder"]),
        ref=config["reference"],
        recal="%s/GATK/recal/${OutputFile}.recal.table" % (config["project-folder"])
    output:
        "%s/BAM/"${OutputFile}-pe.dedub.recal.bam" % (config["project-folder"])        
    log:
        "%s/logs/GATK/PrintReads_{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/GATK/PrintReads_{samples}.benchmark.tsv" % (config["project-folder"])
    params:
        threads=config["params"]["gatk"]["threads"]
        known=config["known-variants"]
    singularity: config["singularity"]["1kbulls"]
    shell:"""
        java -Xmx80G -jar GATK.jar –T PrintReads –nct {params.threads} -R {input.ref} -I {input.bam} -BQSR {input.recal} -o {output}
    """
      
rule AnalyzeCovariates:
   """
    Analyze Covariates (PICARD)
    """
    input:
        bam="%s/BAM/"${OutputFile}-pe.dedub.bam" % (config["project-folder"]),
        ref=config["reference"],
        recal="%s/GATK/recal/${OutputFile}.recal.table" % (config["project-folder"])
    output:
        table="%s/GATK/recal/after_recal.table" % (config["project-folder"]),
        pdf="%s/GATK/recal/after_recal.table" % (config["project-folder"])
    log:
        "%s/logs/GATK/AnalyzeCovariates_{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/GATK/AnalyzeCovariates_{samples}.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["1kbulls"]
    shell:"""
        java -Xmx80G -jar $GATK.jar –T AnalyzeCovariates -R {input.ref} -before {input.table} -after {output.table} -plots {output.pdf}
    """
 
An example GATK AnalyzeCovariates command
