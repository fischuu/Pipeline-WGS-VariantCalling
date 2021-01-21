rule Concatenate_lanes:
    """
    Concatenate the demultiplexed fastq lane files (BASH).
    """
    input:
      R1=expand("%s/FASTQ/RAW/{rawsamples}_R1_001.fastq.gz" % (config["project-folder"]), rawsamples=rawsamples),
      R2=expand("%s/FASTQ/RAW/{rawsamples}_R2_001.fastq.gz" % (config["project-folder"]), rawsamples=rawsamples)
    output:
      R1=temp("%s/FASTQ/CONCATENATED/{samples}_R1.fastq.gz" % (config["project-folder"])),
      R2=temp("%s/FASTQ/CONCATENATED/{samples}_R2.fastq.gz" % (config["project-folder"]))
    log:
        "%s/logs/Concatenate/catFastq_{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/Concatenate/{samples}.benchmark.tsv" % (config["project-folder"])
    params:
       infolder="%s/FASTQ/RAW" % (config["project-folder"]),
       outfolder="%s/FASTQ/CONCATENATED" % (config["project-folder"])
    shell:"""
        mkdir -p {params.outfolder}
        cat {params.infolder}/{wildcards.samples}*_R1_001.fastq.gz > {output.R1} 2> {log}
        cat {params.infolder}/{wildcards.samples}*_R2_001.fastq.gz > {output.R2} 2> {log}
  	"""    