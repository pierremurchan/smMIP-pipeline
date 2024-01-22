# smMIP-pipeline

---
## Introduction

Pipeline for processing and analysing single-molecule molecular inversion probes data usin [smMIP-tools](https://github.com/abelson-lab/smMIP-tools).

---

## Usage


### Running the pipeline

The pipeline can be run with default parameters by running the following command:

```
nextflow run main.nf --input ./samplesheet.csv
```

### Input specifications

The pipeline requires a samplesheet in csv format as input. The samplesheet should contain the following columns: `sample,fastq_1,fastq2`

An example samplesheet is given below:

|    sample   |          fastq_1        |         fastq_2         |
|-------------|-------------------------|-------------------------|
| Sample_ID_1 | Sample_ID_1_R1.fastq.gz | Sample_ID_1_R2.fastq.gz |
| Sample_ID_2 | Sample_ID_2_R1.fastq.gz | Sample_ID_2_R2.fastq.gz |
| Sample_ID_3 | Sample_ID_3_R1.fastq.gz | Sample_ID_3_R2.fastq.gz |


A phenotype file is also required as input to describe the configuration of the assay. The phenotype file should be a tsv file with the following columns: `id type    replicate`

An annotate smMIP design file can also be passed to the pipeline. If an annotated smMIP design file is not provided, a smMIP design file generated with [MIPgen](https://shendurelab.github.io/MIPGEN/) can be passed to the pipeline and the pipeline will attempt to annotate the smMIP design file using the [cellbaseR](https://bioconductor.org/packages/release/bioc/html/cellbaseR.html).