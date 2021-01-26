rule Align_data:
    """
    Align the data (bwa).
    """
    input:
        R1="%s/FASTQ/TRIMMED/{rawsamples}_R1.fastq.gz" % (config["project-folder"]),
        R2="%s/FASTQ/TRIMMED/{rawsamples}_R2.fastq.gz" % (config["project-folder"]),
        ref=config["reference"]
    output: 
        temp("%s/SAM/{rawsamples}-pe.sam" % (config["project-folder"]))
    log:
        "%s/logs/Bwa/alignFastq_{rawsamples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/Bwa/{rawsamples}.benchmark.tsv" % (config["project-folder"])
    params:
        rgid = lambda wildcards: list(samplesheet.lane[samplesheet.rawsample == wildcards.rawsamples]),
        rgpl = config["params"]["bwa"]["rgpl"],
        rgsm = lambda wildcards: list(samplesheet.intid[samplesheet.rawsample == wildcards.rawsamples])
    singularity: config["singularity"]["1kbulls"]
    shell:"""
          bwa mem -M -t 12 -R @RG\\tID:{params.rgid}\\tPL:{params.rgpl}\\tSM:{params.rgsm} \
          {input.ref} {input.R1} {input.R2} > {output}
  	"""
  	
rule sort_and_index:
    """
    Sort and index the alignments (samtools).
    """
    input:
        "%s/SAM/{rawsamples}-pe.sam" % (config["project-folder"])
    output:
        "%s/BAM/{rawsamples}-pe.sorted.bam" % (config["project-folder"])
    log:
        "%s/logs/Samtools/sortIndexFastq_{rawsamples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/Samtools/{rawsamples}.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["1kbulls"]
    shell:"""
        samtools sort -o {output} -O BAM {input}
        samtools index {output}
    """
             
#rule merge_sam_files:
#    """
#    Merge the lane-wise BAM-files into sample-wise BAMs
#    """
#    
#java -Xmx80G -jar  /usr/local/picard/2.1.0/picard.jar MergeSamFiles ${BAMlist} O= ${INTERNATIONALID}.sorted.bam VALIDATION_STRINGENCY=LENIENT #ASSUME_SORTED=true MERGE_SEQUENCE_DICTIONARIES=true


rule mark_duplictes:
    """
    Mark the duplicates in the BAM files
    """
    input:
        "%s/BAM/{samples}-pe.sorted.bam" % (config["project-folder"])
    output:
        bam="%s/BAM/{samples}-pe.dedup.bam" % (config["project-folder"]),
        metric="%s/BAM/{samples}-pe.dedup.metrics" % (config["project-folder"])
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
