#!/usr/bin/env nextflow

/*
 * Define the default parameters - will be moved to external config files later
 */

params.verbose = false // Enable for more verbose information
params.outDir = "${PWD}" // Path to output directory


/*
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/



tsvPath = ''
if (params.sample) tsvPath = params.sample

fastqFiles = Channel.empty()
bamFiles = Channel.empty()
indexRead_file = Channel.value()

tsvFile = file(tsvPath)
fastqFiles = extractFastq(tsvFile)
indexRead_file = extractUMIRead(tsvFile)

(patientGenders, fastqFiles) = extractGenders(fastqFiles)

directoryMap = defineDirectoryMap()
referenceMap = defineReferenceMap()

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

startMessage()

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

  script:
  """
  fastqc -q ${fastqFile1} ${fastqFile2}
  """
}

if (params.verbose) fastQCreport = fastQCreport.view {
  "FastQC report:\n\
  Files : [${it[0].fileName}, ${it[1].fileName}]"
}

process TrimReads {
  tag {idPatient + "-" + idRun}

  input:
    set idPatient, status, idSample, idRun, file(fastqFile1), file(fastqFile2) from fastqFiles

  output:
    file "*fastq.gz" into trimmed_reads
    set idPatient, status, idSample, idRun into trim_output

  script:
  """
  java -Xmx6g -jar /AGeNT/SurecallTrimmer.jar -hs -fq1 ${fastqFile1} -fq2 ${fastqFile2}
  """
}

//Map trimmed reads with BWA. Trim 1bp to reduce error rate. Output into two channels.
process MapReads {
  tag {idPatient + "-" + idRun}

  publishDir "${directoryMap.MapReads}/${idRun}", mode: 'link', pattern: '*trimmed.bam*'
  
  input:
    set idPatient, status, idSample, idRun from trim_output
    file reads from trimmed_reads
    set file(genomeFile), file(bwaIndex) from Channel.value([referenceMap.genomeFile, referenceMap.bwaIndex])

  output:
    set idPatient, status, idSample, idRun, file("${idRun}.bam") into mappedBam
    set idPatient, status, idSample, idRun, file("${idRun}.standard.sorted.trimmed.bam"), file("${idRun}.standard.sorted.trimmed.bam.bai")  into trimmed_StandardBAM

  script:
  readGroup = "@RG\\tID:${idRun}\\tPU:${idRun}\\tSM:${idSample}\\tLB:${idSample}\\tPL:illumina"
  // adjust mismatch penalty for tumor samples
  """
  bwa mem -R \"${readGroup}\" -t ${task.cpus} -B 3 -M \
  ${genomeFile} $reads | \
  samtools sort --threads ${task.cpus} -m 4G - > ${idRun}.bam
  bam trimBam ${idRun}.bam ${idRun}.standard.sorted.trimmed.bam 1
  samtools index ${idRun}.standard.sorted.trimmed.bam
  """
}

process AddUMIs {
  tag {idPatient + "-" + idRun}

  publishDir "${directoryMap.AddUMIs}/${idRun}", mode: 'link'

  input:
    set idPatient, status, idSample, idRun, file(bam) from mappedBam
    file indexRead from indexRead_file

  output:
    set idPatient, status, idSample, idRun, file("${idRun}.UMI.sorted.trimmed.bam"), file("${idRun}.UMI.sorted.trimmed.bam.bai") into trimmed_umiBAM

  script:
  """
  java -Xmx10g -jar /AGeNT/LocatIt.jar -U -IB -b '/Users/viklj600/Projects/TargetSeq/48379-1504167684_Amplicons.bed' -o ${idRun}.UMI.bam ${bam} ${indexRead}
  samtools sort -o ${idRun}.UMI.sorted.bam ${idRun}.UMI.bam
  bam trimBam ${idRun}.UMI.sorted.bam ${idRun}.UMI.sorted.trimmed.bam 1
  samtools index ${idRun}.UMI.sorted.trimmed.bam
  """
}

/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

def checkFileExtension(it, extension) {
  // Check file extension
  if (!it.toString().toLowerCase().endsWith(extension.toLowerCase())) exit 1, "File: ${it} has the wrong extension: ${extension} see --help for more information"
}

def checkParamReturnFile(item) {
  params."${item}" = params.genomes[params.genome]."${item}"
  return file(params."${item}")
}

def defineDirectoryMap() {
  return [
    'bamQC'            : "${params.outDir}/Reports/bamQC",
    'bcftoolsStats'    : "${params.outDir}/Reports/BCFToolsStats",
    'fastQC'           : "${params.outDir}/Reports/FastQC",
    'samtoolsStats'    : "${params.outDir}/Reports/SamToolsStats",
    'MapReads'         : "${params.outDir}/BAMfiles",
    'AddUMIs'          : "${params.outDir}/BAMfiles"
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

      checkFileExtension(fastqFile1,".fastq.gz")
      checkFileExtension(fastqFile2,".fastq.gz")

      [idPatient, gender, status, idSample, idRun, fastqFile1, fastqFile2]
    }
}

def extractUMIRead(tsvFile) {
  // Channeling the TSV file containing FASTQ.
  // Format is: "subject gender status sample lane fastq1 fastq2"
  Channel
    .from(tsvFile.readLines())
    .map{line ->
      def list       = returnTSV(line.split(),8)
      def indexRead = returnFile(list[7])

      checkFileExtension(indexRead,".fastq.gz")

      [indexRead]
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

def minimalInformationMessage() {
  // Minimal information message
  log.info "Command Line: " + workflow.commandLine
  log.info "Profile     : " + workflow.profile
  log.info "Project Dir : " + workflow.projectDir
  log.info "Launch Dir  : " + workflow.launchDir
  log.info "Work Dir    : " + workflow.workDir
  log.info "Out Dir     : " + params.outDir
  log.info "TSV file    : ${tsvFile}"
  log.info "Genome      : " + params.genome
  log.info "Genome_base : " + params.genome_base
  log.info "Reference files used:"
  log.info "  genome      :\n\t" + referenceMap.genomeFile
  log.info "\t" + referenceMap.genomeDict
  log.info "\t" + referenceMap.genomeIndex
  log.info "  bwa indexes :\n\t" + referenceMap.bwaIndex.join(',\n\t')
}

def returnFile(it) {
  // return file if it exists
  if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
  return file(it)
}

def returnTSV(it, number) {
  // return TSV if it has the correct number of items in row
  if (it.size() != number) exit 1, "Malformed row in TSV file: ${it}, see --help for more information"
  return it
}

def returnStatus(it) {
  // Return status if it's correct
  // Status should be only 0 or 1
  // 0 being normal
  // 1 being tumor (or relapse or anything that is not normal...)
  if (!(it in [0, 1])) exit 1, "Status is not recognized in TSV file: ${it}, see --help for more information"
  return it
}

def startMessage() {
  // Display start message
  this.minimalInformationMessage()
}