## SarGet
Hacked version of [Sarek](https://github.com/SciLifeLab/Sarek) for processing of Haloplex HS targeted sequencing data.

Tools are dockerized and will be pulled when needed. For use on Uppmax docker files should be converted to singularity images and stored in a single folder.

##### Steps
1. Trimming with SurecallTrimmer
2. Alignment with BWA
3. UMI mapping with LocatIt
4. Variant calling with Pisces
5. Annotation with snpeff/snpsift
6. Conversion to readable txt file


##### Usage
###### Demultiplexing
Command
```
bcl2fastq --use-bases-mask Y*,I8,Y10,Y*  --minimum-trimmed-read-length 0 --mask-short-adapter-reads 0 --create-fastq-for-index-reads --no-lane-splitting
```

Create sample sheet containing all samples

###### TSV file format
```
idPatient1	 idSample1	 fastqFile1	 fastqFile2	UMI-read
idPatient2	 idSample2	 fastqFile1	 fastqFile2	UMI-read
```
idPatient and idSample is redundant

###### Command
```bash
 nextflow run test_main.nf --sample samplesheet.tsv --ampliconFile ampliconfile.bed --regionsFile regionsfile.bed
```


##### TODO
Enable high sens calling from UMI data. Not working:
```bash
dotnet /Pisces/5.2.7.47/Pisces_5.2.7.47/Pisces.dll -g ${piscesGenome} -bam ${bam} -i ${regions} -OutFolder . -MinVF 0.0005 -SSFilter false -MinBQ 65 -MaxVQ 100 -MinDepthFilter 500 -MinVQ 0 -VQFilter 20 -ReportNoCalls True -CallMNVs False -RMxNFilter 5,9,0.35 -MinDepth 5 -threadbychr true -gVCF false
```
Add coverage analysis again
Add %on_target_ratio
Add percentage of duplicates
Add calling output for non-UMI sample?