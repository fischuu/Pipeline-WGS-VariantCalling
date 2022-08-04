rule Extract_DepthOfCoverage:
    """
    Extract the depth of coverage (GATK)
    """
    input:
        "%s/GATK/DepthOfCoverage/{samples}_dedup_recal.coverage_oneChr.sample_summary" % (config["project-folder"])
    output:
        temp("%s/GATK/DepthOfCoverage/{samples}_dedup_recal.coverage_oneChr.sample_summary.2ndline" % (config["project-folder"]))
    shell:"""
        sed '2q;d' {input} > {output}
    """
    
rule Combine_DepthOfCoverage:
    """
    Combine the depth of coverage (GATK)
    """
    input:
        expand("%s/GATK/DepthOfCoverage/{samples}_dedup_recal.coverage_oneChr.sample_summary.2ndline" % (config["project-folder"]), samples=samples)
    output:
        file="%s/GATK/DepthOfCoverage/Coverage.samples_summary" % (config["project-folder"]),
        md5="%s/GATK/DepthOfCoverage/Coverage.samples_summary.md5" % (config["project-folder"])
    shell:"""
        cat {input} > {output.file}
        md5sum {output.file} > {output.md5}
    """
    
rule R_finalReport:
    """
    Create the final report (R).
    """
    input:
        "%s/GATK/DepthOfCoverage/Coverage.samples_summary" % (config["project-folder"]),
        "%s/GATK/Cohort.g.vcf.gz" % (config["project-folder"]),
        script=config["report-script"]
    output:
        "%s/finalReport.html" % (config["project-folder"])
    log:
        "%s/logs/R/finalReport.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/R/finalReport.benchmark.tsv" % (config["project-folder"])
    singularity:
        config["singularity"]["r-gbs"]
    params:
       projFolder=config["project-folder"],
       pipeFolder=config["pipeline-folder"],
       refGenome=config["reference"],
       samplesheet=config["samplesheet"],
       variants=config["known-variants"]
    shell:"""
       R -e "projFolder <- '{params.projFolder}'; \
             pipelineFolder <- '{params.pipeFolder}'; \
             refGenome.file <- '{params.refGenome}'; \
             samplesheet.file <- '{params.samplesheet}'; \
             variants.file <- '{params.variants}'; \
             snakemake <- TRUE;\
             rmarkdown::render('{input.script}',output_file='{output}')" &> {log}
    """