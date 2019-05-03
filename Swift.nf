#!/usr/bin/env nextflow

/*
 * Define the default parameters - will be moved to external config files later
 */

params.verbose = false // Enable for more verbose information
params.outDir = "${PWD}" // Path to output directory
params.regionsFile = ''
params.masterFile = ''

/*
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

tsvPath = ''
if (params.sample) tsvPath = params.sample

regionsFile = ''
masterFile = ''
masterFile = file(params.masterFile)
regionsFile = file(params.regionsFile)

fastqFiles = Channel.empty()
bamFiles = Channel.empty()

tsvFile = file(tsvPath)
fastqFiles = extractFastq(tsvFile)

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

process TrimReads {
  tag {idPatient}

  input:
    set idPatient, idSample, file(fastqFile1), file(fastqFile2) from fastqFiles

  output:
    file "*fastq.gz" into trimmed_reads
    set idPatient, idSample into trim_output

  script:
  """
    /bbmap/bbduk.sh -Xmx8g in1=${fastqFile1} in2=${fastqFile2} out1=${idSample}_R1_adapter_trimmed.fastq.gz out2=${idSample}_R2_adapter_trimmed.fastq.gz ref=/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
  """
}

//Map trimmed reads with BWA. Output into two channels.
process MapReads {
  tag {idPatient}

  publishDir "${directoryMap.MapReads}/${idSample}", mode: 'link', pattern: '*trimmed.bam*'
  
  input:
    set idPatient, idSample from trim_output
    file reads from trimmed_reads
    file(master) from Channel.value(masterFile)
    set file(genomeFile), file(bwaIndex) from Channel.value([referenceMap.genomeFile, referenceMap.bwaIndex])

  output:
    set idPatient, idSample, file("${idSample}.sorted.trimmed.bam"), file("${idSample}.sorted.trimmed.bam.bai")  into trimmed_BAM

  script:
  readGroup = "@RG\\tID:${idSample}\\tPU:${idSample}\\tSM:${idSample}\\tLB:${idSample}\\tPL:illumina"
  // adjust mismatch penalty for tumor samples
  """
  # Align with BWA
  bwa mem -R \"${readGroup}\" -t ${task.cpus} -B 3 -M -t 32 ${genomeFile} $reads > ${idSample}.sam
  
  # primerclip
  /primerclip/primerclip ${master} ${idSample}.sam ${idSample}_primer_trimmed.sam

  # sort and index
  samtools sort --threads ${task.cpus} -m 4G ${idSample}_primer_trimmed.sam > ${idSample}.sorted.trimmed.bam
  samtools index ${idSample}.sorted.trimmed.bam

  """
}

(trimmed_BAM, trimmed_BAM_forpindel) = trimmed_BAM.into(2)

process Pindel {
  tag {idSample}

  publishDir "${directoryMap.VariantCallingPindel}", mode: 'link'

  input:
    set idPatient, idSample, file(bam), file(bai) from trimmed_BAM_forpindel
    file(genomeFile) from Channel.fromPath(referenceMap.genomeFile)
    file(genomeIndex) from Channel.fromPath(referenceMap.genomeIndex)
    file(regions) from Channel.value(regionsFile)

  output:
    set idPatient, idSample, file("*.vcf") into pindelvariants

  """
  echo '${bam}\t225\t${idSample}' > ${idSample}.txt

  /pindel/pindel -i ${idSample}.txt -f ${genomeFile} -o ${idSample} -c 13 -B 60 -x 2 -T 4

  /pindel/pindel2vcf -p ${idSample}_SI -r ${genomeFile} -R human_g1k_v37_decoy -d 20101123 -e 35 -mc 10 -G -is 5 -c 13
  /pindel/pindel2vcf -p ${idSample}_D -r ${genomeFile} -R human_g1k_v37_decoy -d 20101123 -e 35 -mc 10 -G -is 5 -c 13

  """
}

process VariantCalling {
  tag {idSample}

  publishDir "${directoryMap.VariantCallingPisces}", mode: 'link'
  
  input:
    set idPatient, idSample, file(bam), file(bai) from trimmed_BAM
    file piscesGenome from Channel.value(referenceMap.PiscesReference)
    file regionsFile

  output:
    set idPatient, idSample, file("${idSample}.sorted.trimmed.vcf") into variants
  script:
  """
  dotnet /Pisces/5.2.10.49/Pisces_5.2.10.49/Pisces.dll -CallMNVs false -g $piscesGenome -bam $bam -OutFolder . -gVCF false -i $regionsFile -RMxNFilter 5,9,0.35 -MinDepth 40 --minvq 20 -MinVF 0.005
  """
}

process RunVEP {
  tag {"${vcf}"}

    publishDir "${directoryMap.vep}", mode: 'link'

  input:
    set idPatient, idSample, file(vcf) from variants

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
    'bamQC'                 : "${params.outDir}/Reports/bamQC",
    'bcftoolsStats'         : "${params.outDir}/Reports/BCFToolsStats",
    'fastQC'                : "${params.outDir}/Reports/FastQC",
    'samtoolsStats'         : "${params.outDir}/Reports/SamToolsStats",
    'MapReads'              : "${params.outDir}/BAMfiles",
    'VariantCallingPisces'  : "${params.outDir}/VCFFiles/Pisces",
    'VariantCallingPindel'  : "${params.outDir}/VCFFiles/Pindel",
    'vep'                   : "${params.outDir}/Annotation/VEP",
    'snpeffReports'         : "${params.outDir}/Annotation/snpeffreports",
    'snpeff'                : "${params.outDir}/Annotation/snpeff",
    'txtAnnotate'           : "${params.outDir}/Annotation/txtAnnotate",
    'coverage'              : "${params.outDir}/Reports/Coverage"
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
      def list       = returnTSV(line.split(),4)
      def idPatient  = list[0]
      def idSample   = list[1]
      def fastqFile1 = returnFile(list[2])
      def fastqFile2 = returnFile(list[3])

      checkFileExtension(fastqFile1,".fastq.gz")
      checkFileExtension(fastqFile2,".fastq.gz")

      [idPatient, idSample, fastqFile1, fastqFile2]
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