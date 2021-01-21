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
    threads: config["params"]["fastq"]["threads"]
    params:
        outfolder="%s/QC/RAW/" % (config["project-folder"])
    singularity: "docker://singlecellpipeline/fastqc:v0.0.2"
    shell:"""
        mkdir -p {params.outfolder};
        fastqc -t {threads} -o {params.outfolder} --extract {input.R1} &> {log.R1};
        fastqc -t {threads} -o {params.outfolder} --extract {input.R2} &> {log.R2};
    """
  
rule multiqc_quality_control_raw_data:
    """
    Quality control of fastq files, raw data (MULTIQC).
    """
    input:
        R1=expand("%s/QC/RAW/{rawsamples}_R1_001_fastqc.zip" % (config["project-folder"]), rawsamples=rawsamples),
        R2=expand("%s/QC/RAW/{rawsamples}_R2_001_fastqc.zip" % (config["project-folder"]), rawsamples=rawsamples)
    output:
        R1=directory("%s/QC/RAW/multiqc_R1/" % (config["project-folder"])),
        R2=directory("%s/QC/RAW/multiqc_R2/" % (config["project-folder"]))
    log:
        R1="%s/logs/MULTIQC/multiqc_raw_R1.log" % (config["project-folder"]),
        R2="%s/logs/MULTIQC/multiqc_raw_R2.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/MULTIQC/multiqc_raw_R1.benchmark.tsv" % (config["project-folder"])
    params:
       R1="%s/QC/RAW/*_R1_001_fastqc.zip" % (config["project-folder"]),
       R2="%s/QC/RAW/*_R2_001_fastqc.zip" % (config["project-folder"])
    threads: 12
    singularity: "docker://ewels/multiqc:1.9"
    shell:"""
        multiqc -f --interactive -o {output.R1} {params.R1} &> {log.R1};
        multiqc -f --interactive -o {output.R2} {params.R2} &> {log.R2};
    """