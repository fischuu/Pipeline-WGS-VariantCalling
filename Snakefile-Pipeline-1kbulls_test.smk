# vim: set filetype=sh :
import pandas as pd
from snakemake.utils import validate, min_version

import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

##### load config and sample sheets #####

samplesheet = pd.read_table(config["samplesheet"]).set_index("rawsample", drop=False)
rawsamples=list(samplesheet.rawsample)

workdir: config["project-folder"]

wildcard_constraints:
    rawsamples="|".join(rawsamples),

##### run complete pipeline #####

rule all:
    input:
      expand("%s/SAM/{rawsamples}-pe.sam" % (config["project-folder"]), rawsamples=rawsamples)


rule Align_data:
    """
    Align the data (bwa).
    """
    input:
        R1="%s/FASTQ/RAW/{rawsamples}_R1_001.fastq.gz" % (config["project-folder"]),
        R2="%s/FASTQ/RAW/{rawsamples}_R2_001.fastq.gz" % (config["project-folder"]),
        ref=config["reference"]
    output: 
        "%s/SAM/{rawsamples}-pe.sam" % (config["project-folder"])
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
