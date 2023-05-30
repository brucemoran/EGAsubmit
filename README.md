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

  --sampleCsv     [file]      CSV format, headers: sampleID,read1,read2;
                              sampleID is a unique identifier (e.g. sample_1);
                              read1, read2 are full paths to single gzipped fastqs
                              (e.g. /path/to/sample_1/read_1.fastq.gz)

  or

  --sampleCat     [file]      CSV format, headers: sampleID,dir,ext;
                              sampleID is a unique identifier (e.g. sample_1);
                              dir specifies path to directory in which files to be catted
                              (e.g. /path/to/sample_1)
                              ext is the extension to pattern match for read1 fastqs,
                              if a read2 exists add with a ';' and the pattern to match
                              (e.g. _1.fastq.gz, or _1.fastq.gz;_2.fastq.gz)

  --experiment    [str]       Run name, as per 'experiment' in EGA, used to tag output

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
