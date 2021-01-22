rule Trim_data:
    """
    Trim the concatenated data (trimmomatic).
    """
    input:
      R1="%s/FASTQ/CONCATENATED/{samples}_R1.fastq.gz" % (config["project-folder"]),
      R2="%s/FASTQ/CONCATENATED/{samples}_R2.fastq.gz" % (config["project-folder"])
    output:
      R1="%s/FASTQ/TRIMMED/{samples}_R1.fastq.gz" % (config["project-folder"]),
      R1S="%s/FASTQ/TRIMMED/{samples}_R1-singleton.fastq.gz" % (config["project-folder"]),
      R2="%s/FASTQ/TRIMMED/{samples}_R2.fastq.gz" % (config["project-folder"]),
      R2S="%s/FASTQ/TRIMMED/{samples}_R2-singleton.fastq.gz" % (config["project-folder"]),
      S="%s/FASTQ/TRIMMED/{samples}.summary" % (config["project-folder"])
    log:
        "%s/logs/Trimmomatic/trimFastq_{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/Trimmomatic/{samples}.benchmark.tsv" % (config["project-folder"])
    params:
      threads=config["params"]["trimmomatic"]["threads"],
      faadapter=config["fa-adapter"],
      leading=config["params"]["trimmomatic"]["leading"], 
      trailing=config["params"]["trimmomatic"]["trailing"],  
      slidingwindow=config["params"]["trimmomatic"]["slidingwindow"],
      avgqual=config["params"]["trimmomatic"]["avgqual"],
      minlen=config["params"]["trimmomatic"]["minlen"]
    singularity: config["singularity"]["1kbulls"]
    shell:"""
      java -jar /Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads {params.threads} -summary {output.S} \
      {input.R1} {input.R2} {output.R1} {output.R1S} {output.R2} {output.R2S} \
      MINLEN:{params.minlen} ILLUMINACLIP:{params.faadapter}:2:30:3:1:true LEADING:{params.leading} \
      TRAILING:{params.trailing} SLIDINGWINDOW:{params.slidingwindows} AVGQUAL:{params.avgqual} MINLEN:{params.minlen}
  	"""    