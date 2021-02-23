# EGAsubmit: helper Nextflow for [EGA](https://ega-archive.org/) submission

## Remit
This repo is an attempt to make submission less laborious. It does not automate every aspect of EGA submission. It is designed for fastq, maybe we will include CRAM or other formats but for now it allows single and paired fastq formats.

### Overview
The hardest part of EGA submission (to our novice eyes) is EGAcryptor and linking files and samples using CSVs. This seems like a task that could be suited to Nextflow. We therefore take input as a CSV indicating sampleID, /path/to/read1.fastq.gz, /path/to/read2.fastq.gz and use these as input to EGAcryptor and EGA-formatted CSV-generating processes. These will then be output to the specified directories, leaving you a relatively straightforward clicky job in [the submitter portal](https://ega-archive.org/submitter-portal/)

### Usage:
