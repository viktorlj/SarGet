## SarGet
Hacked version of [Sarek](https://github.com/SciLifeLab/Sarek) for processing of Haloplex HS targeted sequencing data.

All tools in docker viktorlj/targetseq

##### Steps
1. Trimming with SurecallTrimmer
2. Alignment with BWA
3. UMI mapping with LocatIt
4. Variant calling with Pisces
5. Annotation with snpeff/snpsift

##### Usage
```bash
 nextflow run test_main.nf --sample samplesheet.tsv --ampliconFile ampliconfile.bed --regionsFile regionsfile.bed
```

##### TSV file format
```
idPatient	 idSample	 fastqFile1	 fastqFile2	UMI-read
```

##### TODO

Add Uppmax functions / paths

Enable high sens calling from UMI data. Not working:
```bash
dotnet /Pisces/5.2.7.47/Pisces_5.2.7.47/Pisces.dll -g ${piscesGenome} -bam ${bam} -i ${regions} -OutFolder . -MinVF 0.0005 -SSFilter false -MinBQ 65 -MaxVQ 100 -MinDepthFilter 500 -MinVQ 0 -VQFilter 20 -ReportNoCalls True -CallMNVs False -RMxNFilter 5,9,0.35 -MinDepth 5 -threadbychr true -gVCF false
```

Add calling output for non-UMI sample?