# Make Snakemake available (we have it installed in the module bioconda/3)
module load bioconda/3

# Create the rulegraph
snakemake -s ~/git/Pipeline-WGS-VariantCalling/Snakefile-Pipeline-1kbulls.smk \
          --configfile /scratch/project_2002561/MastitisDNA/Pipeline-1kbulls_config.yaml \
          --rulegraph | dot -T png > /scratch/project_2002561/MastitisDNA/workflow.png

snakemake -s ~/git/Pipeline-WGS-VariantCalling/Snakefile-Pipeline-1kbulls.smk \
          -j 300 \
          --use-singularity \
          --singularity-args "-B /scratch:/scratch" \
          --configfile /scratch/project_2002561/MastitisDNA/Pipeline-1kbulls_config.yaml \
          --latency-wait 60 \
          --cluster-config ~/git/Pipeline-WGS-VariantCalling/Pipeline-1kbulls_server_config.yaml \
          --cluster "sbatch -t {cluster.time} --account={cluster.account} --gres=nvme:{cluster.nvme} --job-name={cluster.job-name} --tasks-per-node={cluster.ntasks} --cpus-per-task={cluster.cpus-per-task} --mem-per-cpu={cluster.mem-per-cpu} -p {cluster.partition} -D {cluster.working-directory}" \
          $@
