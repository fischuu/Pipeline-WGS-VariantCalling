# Pipeline-1kbulls
The variant calling pipeline as suggested from the 1k Bulls consortium.

## Installation

To clone into the pipeline type

```
git clone git@github.com:fischuu/Pipeline-1kbulls.git
```

To run the pipeline, it is required to have at least Snakemake version 6.x installed. This
can be done e.g. inside a Conda environment, see createSnakemakeENV.txt

## Usage
The pipeline is designed to start with lane-wise fastq files that are concatenated in the first step.
If this is not required, it can also be started directly from the _samples_ object instead of the _rawsamples_.

The raw input files are expected to follow this naming scheme:

{rawsamples}_R1_001.fastq.gz
and
{rawsamples}_R2_001.fastq.gz

meaning the the {rawsamples} contains all possible files after multiplexing.

For example, this in case we have those files in our folder `FASTQ/RAW`:

```
Sample1_S12_L001_R1_001.fastq.gz
Sample1_S12_L002_R1_001.fastq.gz
Sample1_S12_L001_R2_001.fastq.gz
Sample1_S12_L002_R2_001.fastq.gz
Sample2_S8_L001_R1_001.fastq.gz
Sample2_S8_L002_R1_001.fastq.gz
Sample2_S8_L001_R2_001.fastq.gz
Sample2_S8_L002_R2_001.fastq.gz
```

That means, the file `rawsamples`should contain this:

```
Sample1_S12_L001
Sample1_S12_L002
Sample2_S8_L001
Sample2_S8_L002
```

And the files `samples` should contain this:

```
Sample1_S12
Sample2_S8
```

WARNING! Please ensure that the names are not true subsets from each other! In that case, the concatenating will go wrong!

# Applyable rules

Each step has its own output rule. That means, in addion to the existing rule "all" (the default), you can create also the output files how they are after each step like this

```

```

## Preparations:

Actually, rawsamples and samples are not needed anymore, only sampleSheet.tsv is needed!!!

### Create rawsamples files

First, create the `rawsamples`file. For that, change into the project folder
```
cd $PROJECTFOLDER
```

where `$PROJECTFOLDER` is then the full path to the project folder.

There, you can run (if your FASTQ files are located in `FASTQ/RAW`)

```
find FASTQ/RAW/ -name '*R1_001.fastq.gz' | xargs -n1 basename | sed 's/_R1_001.fastq.gz//g' > rawsamples
```

### Create samples files

find FASTQ/RAW/ -name '*R1_001.fastq.gz' | xargs -n1 basename | cut -d '_' -f1 | sort | uniq > samples

This needs maybe some adjustments in the -f1 part (e.g. f1-2) to meet your file names.

Again, take care that these names are not a subset from each other!!

### Create sampleSheet.tsv

```
cp rawsamples sampleSheet.tsv
```

and then the first column should get a header "rawsample"
### Configuration
Then, the run configuration needs to be adjusted. For that, the recommended way is to copy the file Pipeline-1kbulls_config.yaml from the pipeline folder into the project folder and to adjust there the corresponding paths.

### Start script
Finally, copy the file `run_WGS-Pipeline.sh` from the pipeline folder also to your project folder and adjust the paths.

Then you can start the pipeline just by typing

```
bash run_WGS-Pipeline.sh
```
# Current status
Under development and possible not yet stable for production

## TODO
* Optimizing the download of variants and reference file