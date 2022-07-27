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
        rgid = get_lane,
        rgpl = config["params"]["bwa"]["rgpl"],
        rgsm = get_sample,
        threads = config["params"]["bwa"]["threads"]
    singularity: config["singularity"]["wgs"]
    shell:"""
          echo "RGID:" {params.rgid}
          echo "RGSM:" {params.rgsm}
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
        bam=temp("%s/BAM/{rawsamples}-pe.sorted.bam" % (config["project-folder"])),
        bai=temp("%s/BAM/{rawsamples}-pe.sorted.bam.bai" % (config["project-folder"]))
    log:
        "%s/logs/Samtools/sortIndexFastq_{rawsamples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/Samtools/{rawsamples}.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["wgs"]
    shell:"""
        samtools sort -o {output.bam} -O BAM {input} &> {log}
        samtools index {output.bam}
    """

rule mark_duplicates:
    """
    Mark the duplicates in the BAM files
    """
    input:
        "%s/BAM/{rawsamples}-pe.sorted.bam" % (config["project-folder"])
    output:
        bam=temp("%s/BAM/{rawsamples}-pe.dedup.bam" % (config["project-folder"])),
        metric="%s/BAM/metrics/{rawsamples}-pe.dedup.metrics" % (config["project-folder"])
    log:
        "%s/logs/Picard/Deduplicate_{rawsamples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/Picard/Deplicate_{rawsamples}.benchmark.tsv" % (config["project-folder"])
    params:
        dist=config["params"]["picard"]["distance"]
    singularity: config["singularity"]["wgs"]
    shell:"""
        java -Xmx80G -jar /picard.jar MarkDuplicates I={input} O={output.bam} M={output.metric} \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE={params.dist} CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT &> {log}
    """

rule merge_bam_files:
    """
    Merge the lane-wise BAM-files into sample-wise BAMs
    """
    input:
        get_duplicated_marked_bams
    output:
        bam=temp("%s/BAM/{samples}.sorted.dedup.bam" % (config["project-folder"])),
        bai=temp("%s/BAM/{samples}.sorted.dedup.bam.bai" % (config["project-folder"]))
    log:
        "%s/logs/Picard/merge_{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/Picard/merge_{samples}.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["wgs"]
    params:
        input = lambda wildcards, input: ' '.join('--input ' + v for v in input)
    shell:"""
       java -Xmx80G -jar /picard.jar MergeSamFiles {params.input} O= {output.bam} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true MERGE_SEQUENCE_DICTIONARIES=true &> {log}
       
       samtools index {output.bam}
    """