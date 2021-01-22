rule BaseRecalibration:
   """
    Perform the base recalibration (PICARD)
    """
    input:
        bam="%s/BAM/"${OutputFile}-pe.sorted.bam" % (config["project-folder"]),
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
