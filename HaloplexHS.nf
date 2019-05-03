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
    set idPatient, idSample, file(fastqFile1), file(fastqFile2), file(fastqFile3) from fastqFilesforFastQC

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
    set idPatient, idSample, file(fastqFile1), file(fastqFile2), file(fastqFile3) from fastqFiles

  output:
    file "*fastq.gz" into trimmed_reads
    set idPatient, idSample, file(fastqFile3) into trim_output

  script:
  """
  java -Xmx6g -jar /AGeNT/SurecallTrimmer.jar -hs -fq1 ${fastqFile1} -fq2 ${fastqFile2}
  """
}

//Map trimmed reads with BWA. Trim 1bp to reduce error rate. Output into two channels.
process MapReads {
  tag {idPatient}

  publishDir "${directoryMap.MapReads}", mode: 'link', pattern: '*.bam*'
  
  input:
    set idPatient, idSample, file(fastqFile3) from trim_output
    file reads from trimmed_reads
    set file(genomeFile), file(bwaIndex) from Channel.value([referenceMap.genomeFile, referenceMap.bwaIndex])

  output:
    set idPatient, idSample, file("${idSample}.standard.bam"), file("${idSample}.standard.bam.bai")  into mappedBam

  script:
  readGroup = "@RG\\tID:${idSample}\\tPU:${idSample}\\tSM:${idSample}\\tLB:${idSample}\\tPL:illumina"
  // adjust mismatch penalty for tumor samples
  """
  bwa mem -R \"${readGroup}\" -t ${task.cpus} -B 4 -A 1.0 -w 100 -k 19 -M \
  ${genomeFile} $reads | \
  samtools sort --threads ${task.cpus} -m 4G - > ${idSample}.standard.bam
  samtools index ${idSample}.standard.bam
  """

}

process AddUMIs {
  tag {idPatient}

  publishDir "${directoryMap.AddUMIs}", mode: 'link'

  input:
    set idPatient, idSample, file(bam) from mappedBam
    file(amplicons) from Channel.value(ampliconFile)
    file indexRead from indexRead_file
    file(genomeFile) from Channel.value(referenceMap.genomeFile)
    
  output:
    set idPatient, idSample, file("${idSample}.UMI.sorted.trimmed.bam"), file("${idSample}.UMI.sorted.trimmed.bam.bai") into umiBAM

  script:
  """
  java -Xmx10g -jar /AGeNT/LocatIt.jar -U -IB -m 2 -d 1 -q 25 -b ${amplicons} -o ${idSample}.UMI.bam ${bam} ${indexRead}
  samtools sort -o ${idSample}.UMI.sorted.bam ${idSample}.UMI.bam
  java -jar /fgbio/fgbio-0.8.1.jar ClipBam -r ${genomeFile} -i ${idSample}.UMI.sorted.bam -o ${idSample}.UMI.sorted.trimmed.bam -s --read-one-five-prime=1 --read-one-three-prime=1 --read-two-five-prime=1 --read-two-three-prime=1
  samtools index ${idSample}.UMI.sorted.trimmed.bam
  """
}

// No coverage analysis at the moment
// (trimmed_umiBAM, trimmed_umiBAM_forcoverage) = trimmed_umiBAM.into(2)

process VariantCallingUMI {
  tag {idSample}

  publishDir "${directoryMap.VariantCallingUMI}", mode: 'link'
  
  input:
    set idPatient, idSample, file(bam), file(bai) from umiBAM
    file piscesGenome from Channel.value(referenceMap.PiscesReference)
    file regionsFile

  output:
    set idPatient, idSample, file("${idSample}.UMI.sorted.trimmed.vcf") into variantsUMI
  script:
  """
  dotnet /Pisces/5.2.10.49/Pisces_5.2.10.49/Pisces.dll -CallMNVs false -g $piscesGenome -bam $bam -OutFolder . -gVCF false -i $regionsFile -RMxNFilter 5,9,0.35 -MinDepth 40 --minvq 20 -MinVF 0.005
  """
}

process RunVEP {
  tag {"${vcf}"}

    publishDir "${directoryMap.vep}", mode: 'link'

  input:
    set idPatient, idSample, file(vcf) from variantsUMI

  output:
    file("${vcf.simpleName}_VEP.summary.html") into vepReport
    set idPatient, idSample, file("${vcf.simpleName}_VEP.ann.vcf") into vepVCF

  script:
  genome = params.genome == 'smallGRCh37' ? 'GRCh37' : params.genome
  """
  vep \
  -i ${vcf} \
  -o ${vcf.simpleName}_VEP.ann.vcf \
  --assembly ${genome} \
  --cache \
  --database \
  --everything \
  --fork ${task.cpus} \
  --format vcf \
  --offline \
  --per_gene \
  --stats_file ${vcf.simpleName}_VEP.summary.html \
  --total_length \
  --vcf
  """
}

process siftAddCosmic {
    tag {vcf}

    input:
       set idPatient, idSample, file(vcf) from vepVCF
       set file(cosmic), file(cosmicIndex) from Channel.value([
       referenceMap.cosmic,
       referenceMap.cosmicIndex,
    ])
    
    output:
        set idPatient, idSample, file("${vcf.baseName}.cosmic.ann.vcf") into filteredcosmicvcf

    script:
    """
    java -Xmx4g \
	  -jar \$SNPEFF_HOME/SnpSift.jar \
	  annotate \
	  -info CNT \
    ${cosmic} \
	  ${vcf} \
	  > ${vcf.baseName}.cosmic.ann.vcf
    """

}

process finishVCF {
    tag {vcf}

    publishDir "${directoryMap.txtAnnotate}", mode: 'link'
    
    input:
        set idPatient, idSample, file(vcf) from filteredcosmicvcf

    output:
        file("${idSample}.anno.txt") into finishedFile

    script:
    """
    pisces2pandas.py -i ${vcf} -s ${idSample} -o ${idSample}.anno.txt
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
    'vep'              : "${params.outDir}/Annotation/VEP",
    'snpeffReports'    : "${params.outDir}/Annotation/snpeffreports",
    'snpeff'           : "${params.outDir}/Annotation/snpeff",
    'txtAnnotate'      : "${params.outDir}/Annotation/txtAnnotate",
    'coverage'         : "${params.outDir}/Reports/Coverage"
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
    // Pisces files
    'PiscesReference'  : checkParamReturnFile("PiscesReference"),
    // cosmic VCF with VCF4.1 header
    'cosmic'           : checkParamReturnFile("cosmic"),
    'cosmicIndex'      : checkParamReturnFile("cosmicIndex"),
    // dbNSFP files
    'dbnsfp'           : checkParamReturnFile("dbnsfp"),
    'dbnsfpIndex'      : checkParamReturnFile("dbnsfpIndex")

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
      def fastqFile3 = returnFile(list[4])

      checkFileExtension(fastqFile1,".fastq.gz")
      checkFileExtension(fastqFile2,".fastq.gz")
      checkFileExtension(fastqFile3,".fastq.gz")

      [idPatient, idSample, fastqFile1, fastqFile2, fastqFile3]
    }
}

// Obsolete - remove?
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