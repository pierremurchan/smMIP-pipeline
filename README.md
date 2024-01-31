# smMIP-pipeline

---
## Introduction

Nextflow pipeline for processing and analysing single-molecule molecular inversion probes data using [smMIP-tools](https://github.com/abelson-lab/smMIP-tools).

---

## Installation

Firstly, clone the repository:

```
git clone https://github.com/pierremurchan/smMIP-pipeline
```

Next, install Go and Singuarity using the [official installation instructions](https://docs.sylabs.io/guides/3.0/user-guide/installation.html).

To build the Singularity container (note that this requires root access):

```
singularity build smmip_pipeline.sif singularity/smmip_container.def
```

Finally, install [Nextflow](https://www.nextflow.io/) using the [official installation instructions](https://www.nextflow.io/docs/latest/getstarted.html).

---

## Usage


### Running the pipeline

Generally, the pipeline can be run with default parameters by running:

```
nextflow run main.nf --input ./samplesheet.csv
```

Example data can be found in the `example_data` directory. To run the pipeline on the example data:

```
nextflow run main.nf --input ./example_data/sample_sheet.csv --design_file ./example_data/supplemental_files/Target_MIPgen.txt --annotated_design_file ./example_data/supplemental_files/annotated_Target_MIPgen.txt --phenotype ./example_data/supplemental_files/configuration.txt --output_dir ./example_data/output -with-singularity ./smmip_pipeline.sif
```

> [!NOTE]
> To run the pipeline on the example data, you will need to modify the sample_sheet.csv file to point to the correct paths for the BAM files.

### Input specifications

The pipeline requires a samplesheet in csv format as input. The samplesheet can take two forms.

The first is for using FASTQ files as input to the pipeline. In this case, the sample sheet should contain the following columns: `sample,fastq_1,fastq2`. If FASTQ files originate from multiple lanes or runs, separate rows shold be included for each lane/run using the same sample ID.

An example samplesheet is given below:

|    sample   |             fastq_1              |             fastq_2              |
|-------------|----------------------------------|----------------------------------|
| Sample_ID_1 | /path/to/Sample_ID_1_R1.fastq.gz | /path/to/Sample_ID_1_R2.fastq.gz |
| Sample_ID_2 | /path/to/Sample_ID_2_R1.fastq.gz | /path/to/Sample_ID_2_R2.fastq.gz |
| Sample_ID_3 | /path/to/Sample_ID_3_R1.fastq.gz | /path/to/Sample_ID_3_R2.fastq.gz |


The second is for using BAM files as input to the pipeline. In this case, the sample sheet should contain the following columns: `sample,bam`

An example samplesheet is given below:

|    sample   |              bam             |
|-------------|------------------------------|
| Sample_ID_1 |   /path/to/Sample_ID_1.bam   |
| Sample_ID_2 |   /path/to/Sample_ID_2.bam   |
| Sample_ID_3 |   /path/to/Sample_ID_3.bam   |

<!-- TO DO: Add support for using FASTQ files per lane and cat FASTQ files directly in the pipline.
-->

A phenotype file is also required as input to describe the configuration of the assay. The phenotype file should be a tsv file with the following columns: `id type    replicate`

<!-- TO DO: Modify the phenotype file to be csv format.
-->

An annotate smMIP design file can also be passed to the pipeline. If an annotated smMIP design file is not provided, a smMIP design file generated with [MIPgen](https://shendurelab.github.io/MIPGEN/) can be passed to the pipeline and the pipeline will attempt to annotate the smMIP design file using the [cellbaseR](https://bioconductor.org/packages/release/bioc/html/cellbaseR.html).

### Output specifications

In the output directory, the pipeline will generate the following directories:

```
<outdir>/
├── fastqc/
│   └── ... (FastQC outputs)
|
├── multiqc/
│   └── ... (MultiQC outputs)
|
├── bwamem2/
│   └── ... (BWA-MEM2 aligned BAM files)
|
├── pipeline_info
│   └── ... (Pipeline info)
|
└── smMIP-tools/
    ├── called_mutations.txt
    |
    ├── cleaned_bams/
    │   ├── SAMPLE_ID_1/
    │   |   └── SAMPLE_ID_1_clean.bam
    |   |
    |   ├── SAMPLE_ID_2/
    |   |    └── SAMPLE_ID_2_clean.bam 
    |   └── ...
    |
    └── pileups/
        ├── SAMPLE_ID_1_raw_pileup.txt
        ├── SAMPLE_ID_1_sscs_pileup.txt
        ├── SAMPLE_ID_2_raw_pileup.txt
        ├── SAMPLE_ID_2_sscs_pileup.txt
        └── ...
```



       

