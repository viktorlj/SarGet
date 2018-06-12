#!/usr/bin/env nextflow

/*
 * Define the default parameters - will be moved to external config files later
 */

params.verbose = false // Enable for more verbose information
params.outDir = "${PWD}" // Path to output directory
params.ampliconFile = ''
params.regionsFile = ''

/*
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

tsvPath = ''
if (params.sample) tsvPath = params.sample

ampliconFile = ''
regionsFile = ''
ampliconFile = file(params.ampliconFile)
regionsFile = file(params.regionsFile)

fastqFiles = Channel.empty()
bamFiles = Channel.empty()
indexRead_file = Channel.value()

tsvFile = file(tsvPath)
fastqFiles = extractFastq(tsvFile)
indexRead_file = extractUMIRead(tsvFile)

directoryMap = defineDirectoryMap()
referenceMap = defineReferenceMap()

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

startMessage()

(fastqFiles, fastqFilesforFastQC) = fastqFiles.into(2)

process RunFastQC {
  tag {idPatient}

  publishDir "${directoryMap.fastQC}/${idSample}", mode: 'link'

  input:
    set idPatient, idSample, file(fastqFile1), file(fastqFile2) from fastqFilesforFastQC

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
  tag {idPatient}

  input:
    set idPatient, idSample, file(fastqFile1), file(fastqFile2) from fastqFiles

  output:
    file "*fastq.gz" into trimmed_reads
    set idPatient, idSample into trim_output

  script:
  """
  java -Xmx6g -jar /AGeNT/SurecallTrimmer.jar -hs -fq1 ${fastqFile1} -fq2 ${fastqFile2}
  """
}

//Map trimmed reads with BWA. Trim 1bp to reduce error rate. Output into two channels.
process MapReads {
  tag {idPatient}

  publishDir "${directoryMap.MapReads}/${idSample}", mode: 'link', pattern: '*trimmed.bam*'
  
  input:
    set idPatient, idSample from trim_output
    file reads from trimmed_reads
    set file(genomeFile), file(bwaIndex) from Channel.value([referenceMap.genomeFile, referenceMap.bwaIndex])

  output:
    set idPatient, idSample, file("${idSample}.bam") into mappedBam
    set idPatient, idSample, file("${idSample}.standard.sorted.trimmed.bam"), file("${idSample}.standard.sorted.trimmed.bam.bai")  into trimmed_StandardBAM

  script:
  readGroup = "@RG\\tID:${idSample}\\tPU:${idSample}\\tSM:${idSample}\\tLB:${idSample}\\tPL:illumina"
  // adjust mismatch penalty for tumor samples
  """
  bwa mem -R \"${readGroup}\" -t ${task.cpus} -B 3 -M \
  ${genomeFile} $reads | \
  samtools sort --threads ${task.cpus} -m 4G - > ${idSample}.bam
  bam trimBam ${idSample}.bam ${idSample}.standard.sorted.trimmed.bam 1
  samtools index ${idSample}.standard.sorted.trimmed.bam
  """
}

process AddUMIs {
  tag {idPatient}

  publishDir "${directoryMap.AddUMIs}/${idSample}", mode: 'link'

  input:
    set idPatient, idSample, file(bam) from mappedBam
    file indexRead from indexRead_file
    file(amplicons) from Channel.value(ampliconFile)

  output:
    set idPatient, idSample, file("${idSample}.UMI.sorted.trimmed.bam"), file("${idSample}.UMI.sorted.trimmed.bam.bai") into trimmed_umiBAM

  script:
  """
  java -Xmx10g -jar /AGeNT/LocatIt.jar -U -IB -b ${amplicons} -o ${idSample}.UMI.bam ${bam} ${indexRead}
  samtools sort -o ${idSample}.UMI.sorted.bam ${idSample}.UMI.bam
  bam trimBam ${idSample}.UMI.sorted.bam ${idSample}.UMI.sorted.trimmed.bam 1
  samtools index ${idSample}.UMI.sorted.trimmed.bam
  """
}


process VariantCallingUMI {
 tag {idPatient}

publishDir "${directoryMap.VariantCallingUMI}/${idSample}", mode: 'link'

input:
  set idPatient, idSample, file(bam), file(bai) from trimmed_umiBAM
  file(piscesGenome) from Channel.fromPath(params.genome_base)
  file(regions) from Channel.value(regionsFile)

output:
  set idPatient, idSample, file("${idSample}.UMI.sorted.trimmed.vcf") into variantsUMI


//   High sens settings suggested by Illumina, no output right now
//   dotnet /Pisces/5.2.7.47/Pisces_5.2.7.47/Pisces.dll -g ${piscesGenome} -bam ${bam} -i ${regions} -OutFolder . -MinVF 0.0005 -SSFilter false -MinBQ 65 -MaxVQ 100 -MinDepthFilter 500 -MinVQ 0 -VQFilter 20 -ReportNoCalls True -CallMNVs False -RMxNFilter 5,9,0.35 -MinDepth 5 -threadbychr true -gVCF false

script:
  """
  dotnet /Pisces/5.2.7.47/Pisces_5.2.7.47/Pisces.dll -CallMNVs false -g ${piscesGenome} -bam ${bam} -i ${regions} -OutFolder . -gVCF false -i ${regions}  -RMxNFilter 5,9,0.35
  """

}

// Defined --cache_version 91 to sync with installed version in docker
process RunVEP {
  tag {idPatient}

  publishDir directoryMap.vep, mode: 'link'

  input:
    set idPatient, idSample, file(vcf) from variantsUMI

  output:
    set file("${idSample}.vep.ann.vcf"), file("${idSample}.vep.summary.html") into vepReport

  
  script:
  """
  egrep "#|PASS" ${vcf} > pass.${vcf}
  vep \
  -i pass.${vcf} \
  -o ${idSample}.vep.ann.vcf \
  --stats_file ${idSample}.vep.summary.html \
  --cache \
  --everything \
  --format vcf \
  --offline \
  --per_gene \
  --fork ${task.cpus} \
  --total_length \
  --cache_version 91 \
  --vcf
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
    'AddUMIs'          : "${params.outDir}/BAMfiles",
    'VariantCallingUMI': "${params.outDir}/VCFFiles",
    'vep'              : "${params.outDir}/Annotation/VEP"
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
  // Format is: "subject sample fastq1 fastq2 UMIread"
  Channel
    .from(tsvFile.readLines())
    .map{line ->
      def list       = returnTSV(line.split(),5)
      def idPatient  = list[0]
      def idSample   = list[1]
      def fastqFile1 = returnFile(list[2])
      def fastqFile2 = returnFile(list[3])

      checkFileExtension(fastqFile1,".fastq.gz")
      checkFileExtension(fastqFile2,".fastq.gz")

      [idPatient, idSample, fastqFile1, fastqFile2]
    }
}

def extractUMIRead(tsvFile) {
  // Channeling the TSV file containing FASTQ.
  // Format is: "subject sample lane fastq1 fastq2 UMIread"
  Channel
    .from(tsvFile.readLines())
    .map{line ->
      def list       = returnTSV(line.split(),5)
      def indexRead = returnFile(list[4])

      checkFileExtension(indexRead,".fastq.gz")

      [indexRead]
    }
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

def startMessage() {
  // Display start message
  this.minimalInformationMessage()
}