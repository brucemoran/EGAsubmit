# EGAsubmit: helper Nextflow for [EGA](https://ega-archive.org/) submission

## Remit
This repo is an attempt to make submission less laborious.

It does not automate every aspect of EGA submission.

It is designed for fastq, maybe we will include CRAM or other formats but for now it allows single and paired fastq formats.

## Overview
The hardest part of EGA submission (to my novice eye) is EGAcryptor and linking files and samples using CSVs.

This seems like a task suited to Nextflow.

We take input as a CSV indicating sampleID, /path/to/read1.fastq.gz, /path/to/read2.fastq.gz and use these as input to EGAcryptor and EGA-formatted CSV-generating processes.

These will then be output to the specified directories, leaving you a relatively straightforward clicky job in [the submitter portal](https://ega-archive.org/submitter-portal/)

## Usage
'''
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

  --email         [str]       Email address to send completion notification

  --test          [str]       Use test data, "single" or "paired"
'''

## Output
The output goes into the `runID` dir created in the dir you launch from.

In `runID/` you will find dirs `CSVs` and `EGAcrypted`, and a shell script `*.sh`.

If you gave your egaBox and egaPass credentials, the shell script can be run as-is and your data should upload to your box. This is the preferred option!

If you gave the egaBox, again you can run the shell script. It will ask for your password, possibly for every file this is untested (if you have a egaBox then give egaPass too?!).

If no egaBox or egaPass is given a template is printed.
