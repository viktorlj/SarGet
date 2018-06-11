/*
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

tsvPath = ''
if (params.sample) tsvPath = params.sample

fastqFiles = Channel.empty()
bamFiles = Channel.empty()

tsvFile = file(tsvPath)
fastqFiles = extractFastq(tsvFile)


/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/




/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

def checkFileExtension(it, extension) {
  // Check file extension
  if (!it.toString().toLowerCase().endsWith(extension.toLowerCase())) exit 1, "File: ${it} has the wrong extension: ${extension} see --help for more information"
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
      def indexFile = returnFile(list[7])

      checkFileExtension(fastqFile1,".fastq.gz")
      checkFileExtension(fastqFile2,".fastq.gz")
      checkFileExtension(indexFile,".fastq.gz")

      [idPatient, gender, status, idSample, idRun, fastqFile1, fastqFile2, indexFile]
    }
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