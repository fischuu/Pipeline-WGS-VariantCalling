rule fastqc_quality_control_raw_data:
    """
    Quality control of fastq files, raw data (FASTQC).
    """
    input:
        R1="%s/FASTQ/RAW/{rawsamples}_R1_001.fastq.gz" % (config["project-folder"]),
        R2="%s/FASTQ/RAW/{rawsamples}_R2_001.fastq.gz" % (config["project-folder"])
    output:
        R1="%s/QC/RAW/{rawsamples}_R1_001_fastqc.zip" % (config["project-folder"]),
        R2="%s/QC/RAW/{rawsamples}_R2_001_fastqc.zip" % (config["project-folder"])
    log:
        R1="%s/logs/FASTQC/fastqc_raw_R1.{rawsamples}.log" % (config["project-folder"]),
        R2="%s/logs/FASTQC/fastqc_raw_R2.{rawsamples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/FASTQC/fastqc_raw_R1.{rawsamples}.benchmark.tsv" % (config["project-folder"])
    threads: config["params"]["fastqc"]["threads"]
    params:
        outfolder="%s/QC/RAW/" % (config["project-folder"])
    singularity: "docker://singlecellpipeline/fastqc:v0.0.2"
    shell:"""
        mkdir -p {params.outfolder};
        fastqc -t {threads} -o {params.outfolder} --extract {input.R1} &> {log.R1};
        fastqc -t {threads} -o {params.outfolder} --extract {input.R2} &> {log.R2};
    """
    
rule fastqc_quality_control_trimmed_data:
    """
    Quality control of fastq files, trimmed data (FASTQC).
    """
    input:
        R1="%s/FASTQ/TRIMMED/{rawsamples}_R1.fastq.gz" % (config["project-folder"]),
        R2="%s/FASTQ/TRIMMED/{rawsamples}_R2.fastq.gz" % (config["project-folder"])
    output:
        R1="%s/QC/TRIMMED/{rawsamples}_R1_fastqc.zip" % (config["project-folder"]),
        R2="%s/QC/TRIMMED/{rawsamples}_R2_fastqc.zip" % (config["project-folder"])
    log:
        R1="%s/logs/FASTQC/fastqc_trimmed_R1.{rawsamples}.log" % (config["project-folder"]),
        R2="%s/logs/FASTQC/fastqc_trimmed_R2.{rawsamples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/FASTQC/fastqc_trimmed_R1.{rawsamples}.benchmark.tsv" % (config["project-folder"])
    threads: config["params"]["fastqc"]["threads"]
    params:
        outfolder="%s/QC/TRIMMED/" % (config["project-folder"])
    singularity: "docker://singlecellpipeline/fastqc:v0.0.2"
    shell:"""
        mkdir -p {params.outfolder};
        fastqc -t {threads} -o {params.outfolder} --extract {input.R1} &> {log.R1};
        fastqc -t {threads} -o {params.outfolder} --extract {input.R2} &> {log.R2};
    """