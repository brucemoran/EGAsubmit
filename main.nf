#!/usr/bin/env nextflow

def helpMessage() {
  log.info"""
  -----------------------------------------------------------------------
                            EGAsubmit Nextflow
  -----------------------------------------------------------------------
  Usage:

    nextflow run brucemoran/EGAsubmit

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

    --email         [str]       Email address to send completion notification

    --test          [str]       Use test data, "single" or "paired"
    """.stripIndent()
}

if (params.help) exit 0, helpMessage()

if(!params.test){
  //Test Mandatory Arguments
  if(params.sampleCsv && params.sampleCat){
    exit 1, "Please include only one of --sampleCsv or --sampleCat, see --help for format"
  }

  if(params.sampleCsv == null && params.sampleCat == null){
    exit 1, "Please include one of --sampleCsv or --sampleCat, see --help for format"
  }

  if(!Channel.from(params.experiment, checkIfExists: true)){
      exit 1, "Please include --experiment <your_experiment>"
  }

  if(!Channel.from(params.fastqType, checkIfExists: true)){
      exit 1, "Please include --fastqType <paired or single>"
  }

  //Global variables based on input
  params.outDir = "${params.experiment}"
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

//if we need to cat a sample together
if(params.sampleCat){
  Channel.fromPath("${params.sampleCat}")
         .splitCsv( header: true )
         .map { row -> [row.sampleID, row.dir, row.ext] }
         .set { samplecating }

  process Samplecat {

    label 'low_mem'
    publishDir "${params.outDir}/samples/${sampleID}/cat", mode: "copy"

    input:
    tuple val(sampleID), val(dir), val(ext) from samplecating

    output:
    tuple val(sampleID), file(read1), file(read2) into egac

    script:
    rd1ext = "${ext}".split(';')[0]
    rd2ext = "${ext}".split(';')[1]
    read1 = "${sampleID}.R1.fastq.gz"
    read2 = "${sampleID}.R2.fastq.gz"
    """
    #! bash
    cat \$(find ${dir} | grep ${rd1ext} | sort) > ${read1}
    cat \$(find ${dir} | grep ${rd2ext} | sort) > ${read2}
    """
  }
}

// 0.00: Input using sample.csv, EGAcryptor
if(params.sampleCsv){
  Channel.fromPath("${params.sampleCsv}")
         .splitCsv( header: true )
         .map { row -> [row.sampleID, file(row.read1), file(row.read2)] }
         .set { egac }
}

if(params.fastqType == "paired"){
  process Egacrypt_pe {

    label 'low_mem'
    debug true
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

    echo "${sampleID},${read1}.gpg,${read1}.gpg.md5,${read1}.md5,${read2}.gpg,${read2}.gpg.md5,${read2}.md5" > ${sampleID}.lnk.csv
    """
  }

  //collect outputs from above together as need to make single CSV with all
  process Link_csv_pe {

    label 'low_mem'
    publishDir path: "${params.outDir}/CSVs", mode: "copy"

    input:
    file(lnks) from lnk_csv.collect()

    output:
    file("${params.experiment}.link.csv") into send_link

    script:
    """
    echo "'Sample alias','First Fastq File','First Checksum','First Unencrypted checksum','Second Fastq File','Second Checksum','Second Unencrypted checksum'" > ${params.experiment}.link.csv
    ls *lnk.csv | while read LNK; do
      cat \$LNK >> ${params.experiment}.link.csv
    done
    """
  }
}

if(params.fastqType == "single"){

  process Egacrypt_se {

    label 'low_mem'
    publishDir path: "${params.outDir}/EGAcrypted", mode: "copy", pattern: "*[!csv]"

    input:
    tuple val(sampleID), file(read1) from egac

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

    echo "${sampleID},${read1}.gpg,${read1}.gpg.md5,${read1}.md5" > ${sampleID}.lnk.csv
    """
  }

  //collect outputs from above together as need to make single CSV with all
  process Link_csv_se {

    label 'low_mem'
    publishDir path: "${params.outDir}/CSVs", mode: "copy"

    input:
    file(lnks) from lnk_csv.collect()

    output:
    file("${params.experiment}.link.csv") into send_link

    script:
    """
    echo "'Sample alias','Fastq File','Checksum','Unencrypted checksum'" > ${params.experiment}.link.csv
    ls *lnk.csv | while read LNK; do
      cat \$LNK >> ${params.experiment}.link.csv
    done
    """
  }
}

process Regs_csv {

  label 'low_mem'
  publishDir path: "${params.outDir}/CSVs", mode: "copy"

  input:
  file(regs) from reg_csv.collect()

  output:
  file("${params.experiment}.regs.csv") into ( sh_script, send_regs )

  script:
  """
  echo "title,alias,description,subjectId,bioSampleId,caseOrControl,gender,organismPart,cellLine,region,phenotype" > ${params.experiment}.regs.csv
  ls *reg.csv | while read REG; do
    cat \$REG >> ${params.experiment}.regs.csv
  done
  """
}

process Sh_upload {

  label 'low_mem'
  executor 'local'
  publishDir path: "${params.outDir}", mode: "copy"

  input:
  file(regs) from sh_script

  output:
  file("${params.experiment}.*.sh") into ( sh_script_out, send_scrp )

  script:
  if( params.egaBox && params.egaPass )
    """
    echo "ASPERA_SCP_PASS=${params.egaPass} ascp -P33001  -O33001 -QT -l300M -L- \$(realpath ../../../${params.outDir})/EGAcrypted/* ${params.egaBox}@fasp.ega.ebi.ac.uk:/." > ${params.experiment}.aspera_upload_pass.sh
    """
  else if( params.egaBox && !params.egaPass )
    """
    echo "ascp -P33001 -O33001 -QT -l300M -L- \$(realpath ../../../${params.outDir})/EGAcrypted/* ${params.egaBox}@fasp.ega.ebi.ac.uk:/." > ${params.experiment}.aspera_upload_nopass.sh
    """
  else
    """
    echo "ascp -P33001 -O33001 -QT -l300M -L- \$(realpath ../../../${params.outDir})/EGAcrypted/* <your_ega_box_ID>@fasp.ega.ebi.ac.uk:/." > ${params.experiment}.aspera_upload_nobox_nopass.sh
    """
}

// 4.19: ZIP for sending on sendmail
send_link
  .mix(send_regs)
  .mix(send_scrp)
  .set { sendmail_zip }

process zipup {

  label 'low_mem'
  publishDir "${params.outDir}/zip", mode: 'copy'

  input:
  file(send_all) from sendmail_zip.collect()

  output:
  file("${params.experiment}.EGAsubmit.zip") into send_zip

  script:
  """
  zip -r ${params.experiment}.EGAsubmit.zip *
  """
}

//Completion e-mail notification
if(params.email){
  workflow.onComplete {
    sleep(1000)
    def subject = """\
      [brucemoran/EGAsubmit] SUCCESS: $params.experiment [$workflow.runName]
      """
      .stripIndent()
    if (!workflow.success) {
        subject = """\
          [brucemoran/EGAsubmit] FAILURE: $params.experiment [$workflow.runName]
          """
          .stripIndent()
    }

    def msg = """\
      Pipeline execution summary
      ---------------------------
      experiment  : ${params.experiment}
      RunName     : ${workflow.runName}
      Completed at: ${workflow.complete}
      Duration    : ${workflow.duration}
      workDir     : ${workflow.workDir}
      exit status : ${workflow.exitStatus}
      """
      .stripIndent()

      def attachments = send_zip.toList().getVal()

      sendMail(to: "${params.email}",
               subject: subject,
               body: msg,
               attach: attachments)
  }
}
