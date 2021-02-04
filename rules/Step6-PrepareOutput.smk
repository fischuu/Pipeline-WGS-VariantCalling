rule Extract_DepthOfCoverage:
    """
    Extract the depth of coverage (GATK)
    """
    input:
        "%s/GATK/DepthOfCoverage/{intid}_dedup_recal.coverage.sample_summary" % (config["project-folder"])
    output:
        temp("%s/GATK/DepthOfCoverage/{intid}_dedup_recal.coverage.sample_summary.2ndline" % (config["project-folder"]))
    shell:"""
        sed '2q;d' {input} > {output}
    """
    
rule Combine_DepthOfCoverage:
    """
    Combine the depth of coverage (GATK)
    """
    input:
        expand("%s/GATK/DepthOfCoverage/{intid}_dedup_recal.coverage.sample_summary.2ndline" % (config["project-folder"]), intid=intid)
    output:
        file="%s/GATK/DepthOfCoverage/Coverage.sample_summary" % (config["project-folder"]),
        md5="%s/GATK/DepthOfCoverage/Coverage.sample_summary.md5" % (config["project-folder"])
    shell:"""
        cat {input} > {output}
        md5sum {output.file} > {output.md5}
    """