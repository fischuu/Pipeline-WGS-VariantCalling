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

#rawsamples = pd.read_table(config["samplesheet"], header=None)[0].tolist()
#intid = pd.read_table(config["samplesheet"], header=None)[1].tolist()
#lane = pd.read_table(config["samplesheet"], header=None)[2].tolist()

samplesheet = pd.read_table(config["samplesheet"]).set_index("rawsample", drop=False)
rawsamples=list(samplesheet.rawsample)


workdir: config["project-folder"]

wildcard_constraints:
    rawsamples="|".join(rawsamples),
#    samples="|".join(samples)

def get_sample_lane(sample):
    """Returns lane for given sample."""
    subset = samplesheet.loc[samplesheet['rawsample'] == sample]
    return list(subset['lane'].unique())
    
def get_sample_intid(sample):
    """Returns available lanes for given sample."""
    subset = samplesheet.loc[samplesheet['rawsample'] == sample]
    return list(subset['intid'].unique())
  
##### run complete pipeline #####

rule all:
    input:
      config["known-variants"],
#      config["reference-index"],
#      expand("%s/FASTQ/TRIMMED/{rawsamples}_R1.fastq.gz" % (config["project-folder"]), rawsamples=rawsamples),
      expand("%s/QC/RAW/{rawsamples}_R1_001_fastqc.zip" % (config["project-folder"]), rawsamples=rawsamples),
      expand("%s/QC/TRIMMED/{rawsamples}_R1_fastqc.zip" % (config["project-folder"]), rawsamples=rawsamples),
      expand("%s/SAM/{rawsamples}-pe.sam" % (config["project-folder"]), rawsamples=rawsamples),
      expand("%s/BAM/{rawsamples}-pe.sorted.bam" % (config["project-folder"]), rawsamples=rawsamples)

### setup report #####
report: "report/workflow.rst"

##### load rules #####
include: "rules/Step1-Preparations.smk"
include: "rules/Step2-Trimming.smk"
include: "rules/Step3-QC.smk"
include: "rules/Step4-Alignment.smk"
include: "rules/Step5-VariantCalling.smk"
