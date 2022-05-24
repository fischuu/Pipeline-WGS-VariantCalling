# vim: set filetype=sh :
import pandas as pd
from snakemake.utils import validate, min_version

import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()
shell.executable("bash")

##### WGS variant calling snakemake pipeline     #####
##### Compliant with the 1k Bull requirements    #####
##### Daniel Fischer (daniel.fischer@luke.fi)    #####
##### Natural Resources Institute Finland (Luke) #####
##### Version: 0.1                               #####

##### set minimum snakemake version #####
min_version("5.24")

##### load config and sample sheets #####

samplesheet = pd.read_table(config["samplesheet"]).set_index("rawsample", drop=False)
rawsamples=list(samplesheet.rawsample)
intid=list(samplesheet.intid)

workdir: config["project-folder"]

##### Complete the input configuration
config["reference-dict"] = config["reference"]+".dict" # Here we need to remove the last file ending from reference still...
config["reference-fai"] =  config["reference"]+".fai"
config["reference-index"] = config["reference"]+".amb"

wildcard_constraints:
    rawsamples="|".join(rawsamples),
    intid="|".join(intid)
  
##### run complete pipeline #####

rule all:
    input:
      config["known-variants"],
      config["reference-index"],
      expand("%s/QC/RAW/{rawsamples}_R1_001_fastqc.zip" % (config["project-folder"]), rawsamples=rawsamples),
      expand("%s/QC/TRIMMED/{rawsamples}_R1_fastqc.zip" % (config["project-folder"]), rawsamples=rawsamples),
      expand("%s/GATK/recal/{intid}_recal_plots.pdf" % (config["project-folder"]), intid=intid),
      expand("%s/GATK/GVCF/{intid}_dedup_recal.g.vcf.gz" % (config["project-folder"]), intid=intid),
      expand("%s/GATK/CallableLoci/{intid}.CallableLoci.bed" % (config["project-folder"]), intid=intid),
      expand("%s/GATK/DepthOfCoverage/{intid}_dedup_recal.coverage.sample_summary" % (config["project-folder"]), intid=intid),
#      "%s/GATK/DepthOfCoverage/Coverage.sample_summary" % (config["project-folder"]),
#      "%s/GATK/Cohort.g.vcf.gz" % (config["project-folder"])

### setup report #####
report: "report/workflow.rst"

##### load rules #####
include: "rules/Step1-Preparations.smk"
include: "rules/Step2-Trimming.smk"
include: "rules/Step3-QC.smk"
include: "rules/Step4-Alignment.smk"
include: "rules/Step5-VariantCalling.smk"
#include: "rules/Step6-PrepareOutput.smk"
