import pandas as pd
from snakemake.utils import validate, min_version

import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()
shell.executable("bash")

##### WGS variant calling snakemake pipeline     #####
##### Daniel Fischer (daniel.fischer@luke.fi)    #####
##### Natural Resources Institute Finland (Luke) #####

##### Version: 0.1
version = "0.1"

##### set minimum snakemake version #####
min_version("6.0")

##### load config and sample sheets #####

samplesheet = pd.read_table(config["samplesheet"]).set_index("rawsample", drop=False)
rawsamples=list(samplesheet.rawsample)
samples=list(set(list(samplesheet.sample_name)))  # Unique list
lane=list(samplesheet.lane)

workdir: config["project-folder"]

##### Complete the input configuration
config["reference-dict"] = os.path.splitext(config["reference"])[0]+".dict"
config["reference-fai"] =  config["reference"]+".fai"
config["reference-index"] = config["reference"]+".amb"
config["report-script"] = config["pipeline-folder"]+"/scripts/workflow-report.Rmd"

wildcard_constraints:
    rawsamples="|".join(rawsamples),
    samples="|".join(samples)

##### input checks #####

# rawdata-folder needs to end with "/", add it if missing:
# project-folder should not end with "/", so remove it

##### input function definitions ######

def get_lane(wildcards):
    output = samplesheet.loc[wildcards.rawsamples][["lane"]]
    return output.tolist()

def get_sample(wildcards):
    output = samplesheet.loc[wildcards.rawsamples][["sample_name"]]
    return output.tolist()

def get_raw_input_fastqs(wildcards):
    reads = samplesheet.loc[wildcards.rawsamples][["read1", "read2"]]
    path = config["rawdata-folder"]
    output = [path + x for x in reads]
    return output

def get_raw_input_read1(wildcards):
    reads = samplesheet.loc[wildcards.rawsamples][["read1"]]
    path = config["rawdata-folder"]
    output = [path + x for x in reads]
    return output

def get_raw_input_read2(wildcards):
    reads = samplesheet.loc[wildcards.rawsamples][["read2"]]
    path = config["rawdata-folder"]
    output = [path + x for x in reads]
    return output

def get_duplicated_marked_bams(wildcards):
    rs = samplesheet.loc[samplesheet["sample_name"] == wildcards.samples]["rawsample"]
    prefix = "-pe.dedup.bam"
    outputPlain = [x + prefix for x in rs]
    path = config["project-folder"] + "/BAM/"
    output = [path + x for x in outputPlain]
    return output
    
def get_duplicated_marked_bams_old(wildcards):
    print("Considering now:"+wildcards.samples)
    samplesheet.set_index("sample_name", inplace=True)
    rs = samplesheet.loc[[wildcards.samples]]["rawsample"]
    prefix = "-pe.dedup.bam"
    outputPlain = [x + prefix for x in rs]
    path = config["project-folder"] + "/BAM/"
    output = [path + x for x in outputPlain]
    return output
    
##### Print some welcoming summary #####

##### Print the welcome screen #####
print("#################################################################################")
print("##### Welcome to the WGS pipeline")
print("##### version: "+version)
print("##### Number of rawsamples: "+str(len(rawsamples)))
print("##### Number of samples: "+str(len(samples)))
    
##### run complete pipeline #####

rule all:
    input:
      config["known-variants"],
      config["reference-index"],
      expand("%s/QC/RAW/{rawsamples}_R1_001_fastqc.zip" % (config["project-folder"]), rawsamples=rawsamples),
      expand("%s/QC/TRIMMED/{rawsamples}_R1_fastqc.zip" % (config["project-folder"]), rawsamples=rawsamples),
      expand("%s/GATK/recal/{samples}_recal_plots.pdf" % (config["project-folder"]), samples=samples),
      expand("%s/GATK/GVCF/{samples}_dedup_recal.g.vcf.gz" % (config["project-folder"]), samples=samples),
      expand("%s/GATK/DepthOfCoverage/{samples}_dedup_recal.coverage.sample_summary" % (config["project-folder"]), samples=samples),
      "%s/GATK/DepthOfCoverage/Coverage.samples_summary" % (config["project-folder"]),
      "%s/GATK/output.vcf.gz" % (config["project-folder"]),
      "%s/GATK/output.vqsr.vcf" % (config["project-folder"]),
      "%s/RESULTS/final_variants.vcf" % (config["project-folder"]),
      "%s/finalReport.html" % (config["project-folder"])

rule preparations:
    input:
      config["known-variants"],
      config["reference-index"],
      config["reference-index"],
      config["reference-dict"],
      config["reference-fai"]

rule trimming:
    input:
      expand("%s/FASTQ/TRIMMED/{rawsamples}.summary" % (config["project-folder"]), rawsamples=rawsamples)

rule qc:
    input:
      expand("%s/QC/RAW/{rawsamples}_R1_001_fastqc.zip" % (config["project-folder"]), rawsamples=rawsamples),
      expand("%s/QC/RAW/{rawsamples}_R2_001_fastqc.zip" % (config["project-folder"]), rawsamples=rawsamples),
      expand("%s/QC/TRIMMED/{rawsamples}_R1_fastqc.zip" % (config["project-folder"]), rawsamples=rawsamples),
      expand("%s/QC/TRIMMED/{rawsamples}_R2_fastqc.zip" % (config["project-folder"]), rawsamples=rawsamples)

rule alignment:
    input:
      expand("%s/BAM/metrics/{rawsamples}-pe.dedup.metrics" % (config["project-folder"]), rawsamples=rawsamples),
      expand("%s/BAM/{samples}.sorted.dedup.bam" % (config["project-folder"]), samples=samples, rawsamples=rawsamples)

rule variant_calling:
    input:
      expand("%s/GATK/recal/{samples}_after_recal.table" % (config["project-folder"]), samples=samples),
      expand("%s/GATK/recal/{samples}_recal_plots.pdf" % (config["project-folder"]), samples=samples),
      expand("%s/GATK/GVCF/{samples}_dedup_recal.g.vcf.gz" % (config["project-folder"]), samples=samples),
      expand("%s/GATK/DepthOfCoverage/{samples}_dedup_recal.coverage.sample_summary" % (config["project-folder"]), samples=samples),
      "%s/GATK/Cohort.g.vcf.gz" % (config["project-folder"]),
      "%s/GATK/output.vcf.gz" % (config["project-folder"]),
      "%s/GATK/output.vqsr.vcf" % (config["project-folder"])

rule prepare_output:
    input:
      "%s/GATK/DepthOfCoverage/Coverage.samples_summary" % (config["project-folder"]),
      "%s/finalReport.html" % (config["project-folder"]),
      "%s/RESULTS/final_variants.vcf" % (config["project-folder"])


### setup report #####
report: "report/workflow.rst"

##### load rules #####
include: "rules/Step1-Preparations.smk"
include: "rules/Step2-Trimming.smk"
include: "rules/Step3-QC.smk"
include: "rules/Step4-Alignment.smk"
include: "rules/Step5-VariantCalling.smk"
include: "rules/Step6-PrepareOutput.smk"
