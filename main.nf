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

  println "We are in test mode!"

  params.outDir = "EGAsubmit_test_${params.test}"

  if(params.test == "paired"){

    process test_pe_setup {

      publishDir path: "${params.outDir}", mode: "copy"
      echo true

      input:
      file(csv) from Channel.fromPath("${projectDir}/data/test/sample.pe.csv")

      output:
      file("sample_pe.test.csv") into start_test

      script:
      """
      sed 's#path#${projectDir}/data/test/#g' ${csv} > sample_pe.test.csv
      echo -e "Run:\\nnextflow run brucemoran/EGAsubmit --sampleCsv ${params.outDir}/sample_pe.test.csv"
      """
    }

  }

  if(params.test == "single"){

    Channel
      .fromPath("${projectDir}/data/test/sample.se.csv")
      .set { se_test_csv}

    process test_se_setup {

      publishDir path: "${params.outDir}", mode: "copy"
      echo true

      input:
      file(csv) from se_test_csv

      output:
      file("sample_se.test.csv") into start_test

      script:
      """
      sed 's#path#${projectDir}/data/test/#g' ${csv} > sample_se.test.csv
      echo -e "Run:\\nnextflow run brucemoran/EGAsubmit --sampleCsv ${params.outDir}/sample_se.test.csv --fastqType single"
      """
    }
  }
}

// 0.00: Input using sample.csv, EGAcryptor
if(params.sampleCsv){
  if(params.fastqType == "paired"){

    Channel.fromPath("${params.sampleCsv}")
           .splitCsv( header: true )
           .map { row -> [row.sampleID, file(row.read1), file(row.read2)] }
           .set { egac }

    process egacrypt_pe {

      label 'low_mem'
      publishDir path: "${params.outDir}/EGAcrypted", mode: "copy", pattern: "*[!csv]"

      input:
      tuple val(sampleID), file(read1), file(read2) from egac

      output:
      tuple val(sampleID), file("${read1}.gpg"), file("${read1}.gpg.md5"), file("${read1}.md5"),
      file("${read2}.gpg"), file("${read2}.gpg.md5"), file("${read2}.md5") into files_out
      file("${sampleID}.reg.csv") into reg_csv
      file("${sampleID}.lnk.csv") into lnk_csv

      script:
      """
      java -jar /usr/local/jar/EGAcryptor.jar \
        -i "${read1},${read2}" \
        -o ./

      echo ",${sampleID},,${sampleID},,NA,unknown,,,," > ${sampleID}.reg.csv

      echo "${sampleID},${read1}.gpg,${read1}.gpg.md5,${read1}.md5,${read2}.gpg, ${read2}.gpg.md5,${read2}.md5" > ${sampleID}.lnk.csv
      """
    }

    //collect outputs from above together as need to make single CSV with all

    process link_csv_pe {

      label 'low_mem'
      publishDir path: "${params.outDir}/CSVs", mode: "copy"

      input:
      file(lnks) from lnk_csv.collect()

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
      publishDir path: "${params.outDir}/EGAcrypted", mode: "copy", pattern: "*[!csv]"

      input:
      tuple val(sampleID), file(read1), file(read2) from egac

      output:
      tuple val(sampleID), file("${read1}.gpg"), file("${read1}.gpg.md5"), file("${read1}.md5") into files_out
      file("${sampleID}.reg.csv") into reg_csv
      file("${sampleID}.lnk.csv") into lnk_csv

      script:
      """
      java -jar /usr/local/jar/EGAcryptor.jar \
        -i "${read1}" \
        -o ./

      echo ",${sampleID},,${sampleID},,NA,unknown,,,," > ${sampleID}.reg.csv

      echo "${sampleID},${read1}.gpg,${read1}.gpg.md5,${read1}.md5 > ${sampleID}.lnk.csv
      """
    }

    //collect outputs from above together as need to make single CSV with all
    process link_csv_se {

      label 'low_mem'
      publishDir path: "${params.outDir}/CSVs", mode: "copy"

      input:
      file(lnks) from lnk_csv.collect()

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
    file(regs) from reg_csv.collect()

    output:
    file("${params.runID}.regs.csv")

    script:
    """
    echo "title,alias,description,subjectId,bioSampleId,caseOrControl,gender,organismPart,cellLine,region,phenotype" > ${params.runID}.regs.csv
    ls *reg.csv | while read REG; do
      cat \$REG >> ${params.runID}.regs.csv
    done
    """
  }

  //Completion e-mail notification
  if(params.email){
    workflow.onComplete {
      sleep(100000)
      def subject = """\
        [brucemoran/EGAsubmit] SUCCESS: $params.runID [$workflow.runName]
        """
        .stripIndent()
      if (!workflow.success) {
          subject = """\
            [brucemoran/EGAsubmit] FAILURE: $params.runID [$workflow.runName]
            """
            .stripIndent()
      }

      def msg = """\
        Pipeline execution summary
        ---------------------------
        RunID       : ${params.runID}
        RunName     : ${workflow.runName}
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()

      sendMail(to: "${params.email}",
               subject: subject,
               body: msg)
    }
  }
}
