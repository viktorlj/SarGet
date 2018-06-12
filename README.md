## SarGet
Hacked version of [Sarek](https://github.com/SciLifeLab/Sarek) for processing of Haloplex HS targeted sequencing data.

##### Steps

1. Trimming with SurecallTrimmer
2. Alignment with BWA
3. UMI mapping with LocatIt
4. Variant calling with Pisces

##### Usage
```bash
 nextflow run test_main.nf --sample samplesheet.tsv --ampliconFile ampliconfile.bed --regionsFile regionsfile.bed
```

##### TSV file format
```
idPatient	 idSample	 fastqFile1	 fastqFile2	UMI-read
```

