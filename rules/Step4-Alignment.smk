rule Align_data:
    """
    Align the data (bwa).
    """
    input:
        R1="%s/FASTQ/TRIMMED/{rawsamples}_R1.fastq.gz" % (config["project-folder"]),
        R2="%s/FASTQ/TRIMMED/{rawsamples}_R2.fastq.gz" % (config["project-folder"]),
        ref=config["reference"],
        index=config["reference-index"]
    output: 
        temp("%s/SAM/{rawsamples}-pe.sam" % (config["project-folder"]))
    log:
        "%s/logs/Bwa/alignFastq_{rawsamples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/Bwa/{rawsamples}.benchmark.tsv" % (config["project-folder"])
    params:
        rgid = lambda wildcards: list(samplesheet.lane[samplesheet.rawsample == wildcards.rawsamples]),
        rgpl = config["params"]["bwa"]["rgpl"],
        rgsm = lambda wildcards: list(samplesheet.intid[samplesheet.rawsample == wildcards.rawsamples]),
        threads = config["params"]["bwa"]["threads"]
    singularity: config["singularity"]["1kbulls"]
    shell:"""
          bwa mem -M -t {params.threads} -R \"@RG\\tID:{params.rgid}\\tPL:{params.rgpl}\\tSM:{params.rgsm}\" \
          {input.ref} {input.R1} {input.R2} > {output} 2> {log}
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

rule mark_duplictes:
    """
    Mark the duplicates in the BAM files
    """
    input:
        "%s/BAM/{rawsamples}-pe.sorted.bam" % (config["project-folder"])
    output:
        bam="%s/BAM/{rawsamples}-pe.dedup.bam" % (config["project-folder"]),
        metric="%s/BAM/{rawsamples}-pe.dedup.metrics" % (config["project-folder"])
    log:
        "%s/logs/Picard/Deduplicate_{rawsamples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/Picard/Deplicate_{rawsamples}.benchmark.tsv" % (config["project-folder"])
    params:
        dist=config["params"]["picard"]["distance"]
    singularity: config["singularity"]["1kbulls"]
    shell:"""
        java -Xmx80G -jar /picard.jar MarkDuplicates I={input} O={output.bam} M={output.metric} \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE={params.dist} CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT &> {log}
    """


# THIS IS REALLY BADLY HARD-CODED BUT FOR NOW A QUICK FIX!!!! USE THE BELOW STARTED FUNCTION LATER!!!

rule merge_bam_files:
    """
    Merge the lane-wise BAM-files into sample-wise BAMs
    """
    input:
        de=["%s/BAM/{intid}_L001-pe.dedup.bam" % (config["project-folder"]),
            "%s/BAM/{intid}_L002-pe.dedup.bam" % (config["project-folder"]),
            "%s/BAM/{intid}_L003-pe.dedup.bam" % (config["project-folder"]),
            "%s/BAM/{intid}_L004-pe.dedup.bam" % (config["project-folder"])],
        fake=expand("%s/BAM/{rawsamples}-pe.sorted.bam" % (config["project-folder"]), rawsamples=rawsamples)
    output:
        "%s/BAM/{intid}.sorted.dedup.bam" % (config["project-folder"])
    log:
        "%s/logs/Picard/merge_{intid}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/Picard/merge_{intid}.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["1kbulls"]
    params: " I=".join( ["%s/BAM/{intid}_L001-pe.dedup.bam" % (config["project-folder"]), \
                          "%s/BAM/{intid}_L002-pe.dedup.bam" % (config["project-folder"]), \
                          "%s/BAM/{intid}_L003-pe.dedup.bam" % (config["project-folder"]), \
                          "%s/BAM/{intid}_L004-pe.dedup.bam" % (config["project-folder"])])
    shell:"""
       java -Xmx80G -jar  /picard.jar MergeSamFiles I={params} O= {output} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true MERGE_SEQUENCE_DICTIONARIES=true &> {log}
    """

