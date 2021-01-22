rule Align_data:
    """
    Align the data (bwa).
    """
    input:
        R1="%s/FASTQ/TRIMMED/{samples}_R1.fastq.gz" % (config["project-folder"]),
        R2="%s/FASTQ/TRIMMED/{samples}_R2.fastq.gz" % (config["project-folder"])
        ref=config["reference]
    output:
        temp("%s/SAM/"${OutputFile}-pe.sam" % (config["project-folder"]))
    log:
        "%s/logs/Bwa/alignFastq_{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/Bwa/{samples}.benchmark.tsv" % (config["project-folder"])
    params:
        rgid="LANE NUMBER FROM {RAWSAMPLE}",
        rgpl=config["params"]["bwa"]["rgpl"],
        rgsm=config["params"]["bwa"]["rgsm"]
    singularity: config["singularity"]["1kbulls"]
    shell:"""
          bwa mem -M -t 12 -R @RG\\tID:{parama.rgid}\\tPL:{params.rgpl}\\tSM:{params.rgsm} \
          {input.ref} {input.R1} {input.R2} > {output}
  	"""
  	
rule sort_and_index:
    """
    Sort and index the alignments (samtools).
    """
    input:
        "%s/SAM/"${OutputFile}-pe.sam" % (config["project-folder"])
    output:
        "%s/BAM/"${OutputFile}-pe.sorted.bam" % (config["project-folder"])
    log:
        "%s/logs/Samtools/sortIndexFastq_{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/Samtools/{samples}.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["1kbulls"]
    shell:"""
        samtools sort -o {output} -O BAM {input}
        samtools index {output}
    """
             
rule mark_duplictes:
    """
    Mark the duplicates in the BAM files
    """
    input:
        "%s/BAM/"${OutputFile}-pe.sorted.bam" % (config["project-folder"])
    output:
        bam="%s/BAM/"${OutputFile}-pe.dedup.bam" % (config["project-folder"]),
        metric="%s/BAM/"${OutputFile}-pe.dedup.metrics" % (config["project-folder"])
    log:
        "%s/logs/Picard/Deduplicate_{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/Picard/Deplicate_{samples}.benchmark.tsv" % (config["project-folder"])
    params:
        dist=config["params"]["picard"]["distance"]
    singularity: config["singularity"]["1kbulls"]
    shell:"""
        java -Xmx80G -jar /usr/local/picard/2.18.2/picard.jar MarkDuplicates I={input} O={output.bam} M={output.metric} \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE={params.dist} CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT
    """
