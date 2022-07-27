rule merge_bam_files_old:
    """
    Merge the lane-wise BAM-files into sample-wise BAMs
    """
    input:
        expand("%s/BAM/{rawsamples}-pe.sorted.bam" % (config["project-folder"]), rawsamples=rawsamples)
    output:
        bam=temp("%s/BAM/{samples}.sorted.dedup.bam" % (config["project-folder"])),
        bai=temp("%s/BAM/{samples}.sorted.dedup.bam.bai" % (config["project-folder"]))
    log:
        "%s/logs/Picard/merge_{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/Picard/merge_{samples}.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["1kbulls"]
    params: " I=".join( ["%s/BAM/{samples}_L001-pe.dedup.bam" % (config["project-folder"]), \
                          "%s/BAM/{samples}_L002-pe.dedup.bam" % (config["project-folder"]), \
                          "%s/BAM/{samples}_L003-pe.dedup.bam" % (config["project-folder"]), \
                          "%s/BAM/{samples}_L004-pe.dedup.bam" % (config["project-folder"])])
    shell:"""
       java -Xmx80G -jar  /picard.jar MergeSamFiles I={params} O= {output.bam} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true MERGE_SEQUENCE_DICTIONARIES=true &> {log}
       
       samtools index {output.bam}
    """

