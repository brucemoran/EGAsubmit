# EGAsubmit: helper Nextflow for [EGA](https://ega-archive.org/) submission

## Remit
This repo is an attempt to make submission less laborious.

It does not automate every aspect of EGA submission.

It is designed for fastq, maybe we will include CRAM or other formats but for now it allows single and paired fastq formats.

## Overview
The hardest part of EGA submission (to my novice eye) is EGAcryptor and linking files and samples using CSVs.

This seems like a task suited to Nextflow.

We take input as a CSV specifying sampleID and fastq file(s) and use these as input to EGAcryptor and EGA-formatted CSV-generating processes.

We also output a shell script to allow upload by users. This can also be run when parameter `upload` is specified.

You will then be left with a relatively straightforward clicky job in [the submitter portal](https://ega-archive.org/submitter-portal/)

## Usage
```
nextflow run brucemoran/EGAsubmit -profile singularity

Mandatory parameters:

  -profile        [str]       Configuration profile (required: singularity)

  --sampleCsv     [file]      CSV format, headers: sampleID,read1 (option:
                              add read2 for "paired")

  --runID         [str]       Name for run, used to tag output (tip: use the
                              same as the 'experiment' name in submitter portal)

  --fastqType     [str]       Either "single" or "paired" (default)

Optional parameters:

  --egaBox        [str]       ega-box ID (e.g. ega_box_12345)

  --egaPass       [str]       ega-box password

  --upload        [bool]      Upload to egaBox? Requires egaBox and egaPass

  --email         [str]       Email address to send completion notification

  --test          [str]       Use test data, "single" or "paired"
```

## Output and Uploading
The output goes into the `runID` dir created in the dir you launch from.

In `runID/` you will find dirs `CSVs` and `EGAcrypted`, and a shell script `*.sh` based on your supplied egaBox, egaPass.

If you specify `upload` you will also trigger uploading of data to the egaBox specified. NB otherwise you can shell into the singularity container and upload from there if `aspera` is not installed on your system.
