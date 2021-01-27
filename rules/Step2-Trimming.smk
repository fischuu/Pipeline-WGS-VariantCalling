rule Trim_data:
    """
    Trim the concatenated data (trimmomatic).
    """
    input:
      R1="%s/FASTQ/RAW/{rawsamples}_R1_001.fastq.gz" % (config["project-folder"]),
      R2="%s/FASTQ/RAW/{rawsamples}_R2_001.fastq.gz" % (config["project-folder"])
    output:
      R1="%s/FASTQ/TRIMMED/{rawsamples}_R1.fastq.gz" % (config["project-folder"]),
      R1S="%s/FASTQ/TRIMMED/{rawsamples}_R1-singleton.fastq.gz" % (config["project-folder"]),
      R2="%s/FASTQ/TRIMMED/{rawsamples}_R2.fastq.gz" % (config["project-folder"]),
      R2S="%s/FASTQ/TRIMMED/{rawsamples}_R2-singleton.fastq.gz" % (config["project-folder"]),
      S="%s/FASTQ/TRIMMED/{rawsamples}.summary" % (config["project-folder"])
    log:
        "%s/logs/Trimmomatic/trimFastq_{rawsamples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/Trimmomatic/{rawsamples}.benchmark.tsv" % (config["project-folder"])
    params:
      threads=config["params"]["trimmomatic"]["threads"],
      clipping=config["params"]["trimmomatic"]["illuminaclip"],
      leading=config["params"]["trimmomatic"]["leading"], 
      trailing=config["params"]["trimmomatic"]["trailing"],  
      slidingwindow=config["params"]["trimmomatic"]["slidingwindow"],
      avgqual=config["params"]["trimmomatic"]["avgqual"],
      minlen=config["params"]["trimmomatic"]["minlen"]
    singularity: config["singularity"]["1kbulls"]
    shell:"""
      java -jar /Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads {params.threads} -summary {output.S} \
      {input.R1} {input.R2} {output.R1} {output.R1S} {output.R2} {output.R2S} \
      MINLEN:{params.minlen} ILLUMINACLIP:{params.clipping} LEADING:{params.leading} \
      TRAILING:{params.trailing} SLIDINGWINDOW:{params.slidingwindow} AVGQUAL:{params.avgqual} MINLEN:{params.minlen}
  	"""    