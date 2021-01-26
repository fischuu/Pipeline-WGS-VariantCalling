# vim: set filetype=sh :
import pandas as pd
from snakemake.utils import validate, min_version

import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

report: "report/workflow.rst"

##### RNASeq-snakemake pipeline #####
##### Daniel Fischer (daniel.fischer@luke.fi)
##### Natural Resources Institute Finland (Luke)
##### Version: 0.1

##### set minimum snakemake version #####
min_version("5.32")

##### load config and sample sheets #####

rawsamples = pd.read_table(config["samplesheet"], header=None)[0].tolist()

workdir: config["project-folder"]

wildcard_constraints:
    rawsamples="|".join(rawsamples),
#    samples="|".join(samples)

##### run complete pipeline #####

rule all:
    input:
      config["known-variants"],
      config["reference-index"],
      expand("%s/FASTQ/TRIMMED/{rawsamples}_R1.fastq.gz" % (config["project-folder"]), rawsamples=rawsamples),
      expand("%s/QC/RAW/{rawsamples}_R1_001_fastqc.zip" % (config["project-folder"]), rawsamples=rawsamples),
      expand("%s/QC/TRIMMED/{rawsamples}_R1_fastqc.zip" % (config["project-folder"]), rawsamples=rawsamples)

### setup report #####
report: "report/workflow.rst"

##### load rules #####
include: "rules/Step1-Preparations.smk"
include: "rules/Step2-Trimming.smk"
include: "rules/Step3-QC.smk"
include: "rules/Step4-Alignment.smk"
include: "rules/Step5-VariantCalling.smk"
