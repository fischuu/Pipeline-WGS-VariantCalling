rule Download_known_variants:
    """
    Download the suggested reference genome
    """
    output: 
        config["known-variants"]
    log: 
        "%s/logs/Bash/VariantsDownload.log" % (config["project-folder"])
    benchmark: 
        "%s/benchmark/Bash/VariantsDownload.benchmark.tsv" % (config["project-folder"])
    params:
        url=config["urls"]["known-variants"]
    shell:"""
        wget {params.url}
        mv `basename {params.url}` {output}
  	"""    
        
rule Download_reference:
    """
    Download the suggested reference genome
    """
    output: 
        config["reference"]
    log:
        "%s/logs/Bash/ReferenceDownload.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/Bash/ReferenceDownload.benchmark.tsv" % (config["project-folder"])
    params:
       url=config["urls"]["reference"]
    shell:"""
        curl {params.url} > {output}
  	"""    
                                 
rule Index_reference:
    """
    Create the bwa index for the reference genome
    """
    input:
       config["reference"]
    output: 
       config["reference-index"]
    log:
        "%s/logs/Bwa/IndexReference.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/Bwa/IndexReference.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["1kbulls"]
    shell:"""
        bwa index {input} &> {log}
  	"""    
                                 
rule Fai_Index_reference:
    """
    Create the fasta index for the reference genome
    """
    input:
       config["reference"]
    output: 
       fai=config["reference-fai"],
       dict=config["reference-dict"]
    log:
        "%s/logs/Samtools/IndexReference.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/Samtools/IndexReference.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["1kbulls"]
    shell:"""
     java -Xmx40G -jar /picard.jar  CreateSequenceDictionary R={input} O={output.dict}
    samtools faidx {input} &> {log}
  	"""    