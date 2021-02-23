#!/usr/bin/env nextflow

def helpMessage() {
  log.info"""
  -----------------------------------------------------------------------
                            EGAsubmit Nextflow
  -----------------------------------------------------------------------
  Usage:

    nextflow run brucemoran/EGAsubmit

  Mandatory arguments:

    -profile        [str]       Configuration profile (required: singularity)

    --sampleCsv     [file]      CSV format, headers: sampleID,read1,read2

    --runID         [str]       Name for run, used to tag output

    --fastqType     [str]       Either "single" or "paired" (default)

    --email         [str]       Email address to send completion notification

    --test          [str]       Use test data, "single" or "paired"
    """.stripIndent()
}

if (params.help) exit 0, helpMessage()

if(!params.test){
  //Test Mandatory Arguments
  if(!Channel.from(params.sampleCsv, checkIfExists: true)){
    exit 1, "Please include --sampleCsv, see --help for format"
  }

  if(!Channel.from(params.runID, checkIfExists: true)){
      exit 1, "Please include --runID <your_runID>"
  }

  if(!Channel.from(params.fastqType, checkIfExists: true)){
      exit 1, "Please include --fastqType <paired or single>"
  }

  //Global variables based on input
  params.outDir = "${params.runID}"
}

//testing
if(params.test){
  if(params.test == "paired"){

    Channel
      .fromPath("${projectDir}/data/test/sample.pe.csv", type: 'file')
      .set { pe_test_csv}

    process test_se_setup {

      input:
      file(csv) from pe_test_csv

      output:
      file("sampleCsv.pe.csv") into start_test

      script:
      """
      sed 's#path#${projectDir}/data/test/#g' ${csv} > sampleCsv.pe.csv
      """
    }

    params.sampleCsv = start_test
    params.fastqType = "paired"

  }

  if(params.test == "single"){
    Channel
      .fromPath("${projectDir}/data/test/sample.se.csv", type: 'file')
      .set { se_test_csv}

    process test_se_setup {

      input:
      file(csv) from se_test_csv

      output:
      file("sampleCsv.se.csv") into start_test

      script:
      """
      sed 's#path#${projectDir}/data/test/#g' ${csv} > sampleCsv.se.csv
      """
    }

    params.sampleCsv = start_test
    params.fastqType = "single"

  }
}

// 0.00: Input using sample.csv, EGAcryptor
if(params.fastqType == "paired"){
  Channel.fromPath("${params.sampleCsv}")
         .splitCsv( header: true )
         .map { row -> [row.sampleID, file(row.read1), file(row.read2)] }
         .set { egac }

  process egacrypt_pe {

    label 'low_mem'
    publishDir path: "${params.outDir}/EGAcrypted", mode: "copy"

    input:
    tuple val(sampleID), file(read1), file(read2) from egac

    output:
    tuple val(sampleID), file("${read1}.gpg"), file("${read1}.gpg.md5"), file("${read1}.md5"),
    file("${read2}.gpg"), file("${read2}.gpg.md5"), file("${read2}.md5") into files_out
    tuple val(sampleID), file("${sampleID}.reg.csv") into reg_csv
    tuple val(sampleID), file("${sampleID}.lnk.csv") into lnk_csv

    script:
    def taskmem = task.memory == null ? "" : "-Xmx" + javaTaskmem("${task.memory}")
    """
    java ${taskmem} -jar /usr/local/jar/EGAcryptor.jar \
      -i "${read1},${read2}"

    echo ",${sampleID},,${sampleID},,NA,unknown,,,," > ${sampleID}.reg.csv

    echo "${sampleID},${read1}.gpg,${read1}.gpg.md5,${read1}.md5,${read2}.gpg, ${read2}.gpg.md5,${read2}.md5" > ${sampleID}.lnk.csv
    """
  }

  //collect outputs from above together as need to make single CSV with all

  process link_csv_pe {

    label 'low_mem'
    publishDir path: "${params.outDir}/CSVs", mode: "copy"

    input:
    tuple val(sampleID), file(lnks) from lnk_csv.collect()

    output:
    file("${params.runID}.link.csv")

    script:
    """
    echo "Sample alias,First Fastq File,First Checksum,First Unencrypted checksum,Second Fastq File,Second Checksum,Second Unencrypted checksum" > ${params.runID}.link.csv
    ls *lnk.csv | while read LNK; do
      cat \$LNK >> ${params.runID}.link.csv
    done
    """
  }
}

if(params.fastqType == "single"){
  Channel.fromPath("${params.sampleCsv}")
         .splitCsv( header: true )
         .map { row -> [row.sampleID, file(row.read1)] }
         .set { egac }

  process egacrypt_se {

    label 'low_mem'
    publishDir path: "${params.outDir}/EGAcrypted", mode: "copy"

    input:
    tuple val(sampleID), file(read1), file(read2) from egac

    output:
    tuple val(sampleID), file("${read1}.gpg"), file("${read1}.gpg.md5"), file("${read1}.md5") into files_out
    tuple val(sampleID), file("${sampleID}.reg.csv") into reg_csv
    tuple val(sampleID), file("${sampleID}.lnk.csv") into lnk_csv

    script:
    def taskmem = task.memory == null ? "" : "-Xmx" + javaTaskmem("${task.memory}")
    """
    java ${taskmem} -jar /usr/local/jar/EGAcryptor.jar \
      -i "${read1},${read2}"

    echo ",${sampleID},,${sampleID},,NA,unknown,,,," > ${sampleID}.reg.csv

    echo "${sampleID},${read1}.gpg,${read1}.gpg.md5,${read1}.md5 > ${sampleID}.lnk.csv
    """
  }

  //collect outputs from above together as need to make single CSV with all
  process link_csv_se {

    label 'low_mem'
    publishDir path: "${params.outDir}/CSVs", mode: "copy"

    input:
    tuple val(sampleID), file(lnks) from lnk_csv.collect()

    output:
    file("${params.runID}.link.csv")

    script:
    """
    echo "Sample alias,Fastq File,Checksum,Unencrypted checksum" > ${params.runID}.link.csv
    ls *lnk.csv | while read LNK; do
      cat \$LNK >> ${params.runID}.link.csv
    done
    """
  }
}

process regs_csv {

  label 'low_mem'
  publishDir path: "${params.outDir}/CSVs", mode: "copy"

  input:
  tuple val(sampleID), file("${sampleID}.reg.csv") from reg_csv

  output:
  file("${params.runID}.regs.csv")

  script:
  """
  echo "title,alias,description,subjectId,bioSampleId,caseOrControl,gender,organismPart,cellLine	region,phenotype" > ${params.runID}.regs.csv
  ls *reg.csv | while read REG; do
    cat \$REG >> ${params.runID}.regs.csv
  done
  """
}
