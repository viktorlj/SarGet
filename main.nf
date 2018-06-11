
/*
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

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

if (params.verbose) bamFiles = bamFiles.view {
  "BAMs to process:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  Files : [${it[3].fileName}, ${it[4].fileName}]"
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
      def list       = returnTSV(line.split(),7)
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