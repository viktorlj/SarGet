
/*
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

tsvPath = params.sample

fastqFiles = Channel.empty()
bamFiles = Channel.empty()

fastqFiles = extractFastq(tsvFile);
(patientGenders, fastqFiles) = extractGenders(fastqFiles)

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

(fastqFiles, fastqFilesforFastQC) = fastqFiles.into(2)

if (params.verbose) fastqFiles = fastqFiles.view {
  "FASTQs to preprocess:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\tRun   : ${it[3]}\n\
  Files : [${it[4].fileName}, ${it[5].fileName}]"
}

process RunFastQC {
  tag {idPatient + "-" + idRun}

  publishDir "${directoryMap.fastQC}/${idRun}", mode: 'link'

  input:
    set idPatient, status, idSample, idRun, file(fastqFile1), file(fastqFile2) from fastqFilesforFastQC

  output:
    file "*_fastqc.{zip,html}" into fastQCreport

  when: step == 'mapping' && !params.noReports

  script:
  """
  fastqc -q ${fastqFile1} ${fastqFile2}
  """
}

if (params.verbose) fastQCreport = fastQCreport.view {
  "FastQC report:\n\
  Files : [${it[0].fileName}, ${it[1].fileName}]"
}

/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

def extractFastq(tsvFile) {
  // Channeling the TSV file containing FASTQ.
  // Format is: "subject gender status sample lane fastq1 fastq2"
  Channel
    .from(tsvFile.readLines())
    .map{line ->
      def list       = returnTSV(line.split(),8)
      def idPatient  = list[0]
      def gender     = list[1]
      def status     = returnStatus(list[2].toInteger())
      def idSample   = list[3]
      def idRun      = list[4]
      def fastqFile1 = returnFile(list[5])
      def fastqFile2 = returnFile(list[6])
      def indexFile = returnFile(list[7])

      checkFileExtension(fastqFile1,".fastq.gz")
      checkFileExtension(fastqFile2,".fastq.gz")
      checkFileExtension(indexFile,".fastq.gz")

      [idPatient, gender, status, idSample, idRun, fastqFile1, fastqFile2, indexFile]
    }
}

def extractGenders(channel) {
  def genders = [:]  // an empty map
  channel = channel.map{ it ->
    def idPatient = it[0]
    def gender = it[1]
    genders[idPatient] = gender
    [idPatient] + it[2..-1]
  }
  [genders, channel]
}

// Lots of work needed to add correct directories 

def defineDirectoryMap() {
  return [
    'Aligned'     : "${params.outDir}/Preprocessing/NonRealigned",
    'nonRecalibrated'  : "${params.outDir}/Preprocessing/NonRecalibrated",
    'recalibrated'     : "${params.outDir}/Preprocessing/Recalibrated",
    'bamQC'            : "${params.outDir}/Reports/bamQC",
    'bcftoolsStats'    : "${params.outDir}/Reports/BCFToolsStats",
    'fastQC'           : "${params.outDir}/Reports/FastQC",
    'markDuplicatesQC' : "${params.outDir}/Reports/MarkDuplicates",
    'samtoolsStats'    : "${params.outDir}/Reports/SamToolsStats"
  ]
}

def defineReferenceMap() {
  if (!(params.genome in params.genomes)) exit 1, "Genome ${params.genome} not found in configuration"
  return [
    // genome reference dictionary
    'genomeDict'       : checkParamReturnFile("genomeDict"),
    // FASTA genome reference
    'genomeFile'       : checkParamReturnFile("genomeFile"),
    // genome .fai file
    'genomeIndex'      : checkParamReturnFile("genomeIndex"),
    // BWA index files
    'bwaIndex'         : checkParamReturnFile("bwaIndex"),
  ]
}

def checkFileExtension(it, extension) {
  // Check file extension
  if (!it.toString().toLowerCase().endsWith(extension.toLowerCase())) exit 1, "File: ${it} has the wrong extension: ${extension} see --help for more information"
}