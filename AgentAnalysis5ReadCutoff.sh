
now=$('date' +"%Y%m%d")

<<'COMM1'
COMM1
#Trim reads with Agent
java -Xmx6g -jar /Users/viklj600/BioData/AGeNT/SurecallTrimmer_v4.0.1.jar -fq1 ../RawData/"$1"_R1_001.fastq.gz -fq2 ../RawData/"$1"_R2_001.fastq.gz -hs -out_loc temp/

#Rename trimmed reads
for i in temp/*_Cut_0.fastq.gz; do mv $i ${i/fastq.gz*/trimmed.fastq.gz}; done;

mv temp/"$1"_R1_001.trimmed.fastq.gz Trimmed/
mv temp/"$1"_R2_001.trimmed.fastq.gz Trimmed/

#Align
/Users/viklj600/BioData/bwa-0.7.12/bwa mem -t 8 -M -R "@RG\tID:"$now"_${1}\tSM:${1}\tPL:illumina" /Volumes/Cambridge/Reference/hg19.fa Trimmed/"$1"_R1_001.trimmed.fastq.gz Trimmed/"$1"_R2_001.trimmed.fastq.gz | samtools view -hb - | samtools sort -@ 6 -o temp/"$1".sorted.bam -
samtools index temp/"$1".sorted.bam

#Add UMIs to bam file
java -Xmx12G -jar /Users/viklj600/BioData/AGeNT/LocatIt_v4.0.1.jar -X temp -t temp -U -IB -b 48379-1504167684_Amplicons.bed -o temp/"$1".UMI.unsorted.bam temp/"$1".sorted.bam ../RawData/"$1"_UMI_001.fastq.gz

samtools sort -o temp/"$1".UMI.sorted.bam temp/"$1".UMI.unsorted.bam
samtools index temp/"$1".UMI.sorted.bam

#Trim 1bp from both bam-files
/Users/viklj600/BioData/bamUtil_1.0.13/bamUtil/bin/bam trimBam temp/"$1".UMI.sorted.bam BamFiles/"$1".UMI.sorted.trimmed.bam 1
samtools index BamFiles/"$1".UMI.sorted.trimmed.bam

/Users/viklj600/BioData/bamUtil_1.0.13/bamUtil/bin/bam trimBam temp/"$1".sorted.bam BamFiles/"$1".standard.sorted.trimmed.bam 1
samtools index BamFiles/"$1".standard.sorted.trimmed.bam

mv temp/"$1".UMI.unsorted.properties  Histograms/

#Variant calling
samtools mpileup -f /Volumes/Cambridge/Reference/hg19.fa -d 100000 -A -B -q 20 -l 48379-1504167684_Regions.bed BamFiles/"$1".UMI.sorted.trimmed.bam | java -Xmx14g -jar /Users/viklj600/BioData/VarScan/VarScan.v2.3.7.jar mpileup2cns --min-var-freq 0.005 --min-coverage 40 --min-avg-qual 20 --variants 1 --min-reads2 5 --output-vcf 1 --strand-filter 0 > VarScanOutUMI/"$1".UMI.variants.vcf
samtools mpileup -f /Volumes/Cambridge/Reference/hg19.fa -d 100000 -A -B -q 20 -l 48379-1504167684_Regions.bed BamFiles/"$1".standard.sorted.trimmed.bam | java -Xmx14g -jar /Users/viklj600/BioData/VarScan/VarScan.v2.3.7.jar mpileup2cns --min-var-freq 0.005 --min-coverage 40 --min-avg-qual 20 --variants 1 --min-reads2 5 --output-vcf 1 --strand-filter 0 > VarScanOutStandard/"$1".standard.variants.vcf

#Annovar
perl /Users/viklj600/BioData/annovar/convert2annovar.pl -format vcf4 VarScanOutUMI/"$1".UMI.variants.vcf --outfile AnnovarUMI/"$1".UMI.annovarinput -include
perl /Users/viklj600/BioData/annovar/convert2annovar.pl -format vcf4 VarScanOutStandard/"$1".standard.variants.vcf --outfile AnnovarStandard/"$1".standard.annovarinput -include


perl /Users/viklj600/BioData/annovar/table_annovar.pl AnnovarUMI/"$1".UMI.annovarinput /Users/viklj600/BioData/annovar/humandb/ -protocol refGene,1000g2014oct_eur,snp138,snp138NonFlagged,cosmic82 -operation g,f,f,f,f -nastring "-" -otherinfo -buildver hg19 -remove -arg '-splicing_threshold 5',,,,
perl /Users/viklj600/BioData/annovar/table_annovar.pl AnnovarStandard/"$1".standard.annovarinput /Users/viklj600/BioData/annovar/humandb/ -protocol refGene,1000g2014oct_eur,snp138,snp138NonFlagged,cosmic82 -operation g,f,f,f,f -nastring "-" -otherinfo -buildver hg19 -remove -arg '-splicing_threshold 5',,,,


#Convert output
python ConvertAnnovarOutputGermline.py -i AnnovarUMI/"$1".UMI.annovarinput.hg19_multianno.txt -s "$1"_UMI -o AnnotatedUMI/"$1".UMI.annotated.txt
python ConvertAnnovarOutputGermline.py -i AnnovarStandard/"$1".standard.annovarinput.hg19_multianno.txt -s "$1"_standard -o AnnotatedStandard/"$1".standard.annotated.txt

#Coverage analysis
sambamba depth region --filter 'mapping_quality > 20' -q 20 -T 10 -T 50 -T 100 -T 500 -T 1000 -L 48379-1504167684_Regions.bed -o SambambaOutUMI/"$1".UMI.Coverage.txt BamFiles/"$1".UMI.sorted.trimmed.bam
sambamba depth region --filter 'mapping_quality > 20' -q 20 -T 10 -T 50 -T 100 -T 500 -T 1000 -L 48379-1504167684_Regions.bed -o SambambaOutStandard/"$1".standard.Coverage.txt BamFiles/"$1".standard.sorted.trimmed.bam


python ParseSamCoverage.py -i SambambaOutUMI/"$1".UMI.Coverage.txt -o CoverageUMI/"$1".UMI.MeanCoverage.txt
python ParseSamCoverage.py -i SambambaOutStandard/"$1".standard.Coverage.txt -o CoverageStandard/"$1".standard.MeanCoverage.txt

<<'COMM2'

COMM2

#rm temp/*

#UNUSED STEP WHEN DEMUTLIPLEXING WITH bcl2fastq
#Match read1 and read2 with demultiplexed index file
#/Users/viklj600/BioData/ngsutils/bin/fastqutils properpairs -z -f IndexReads/"$1"_I1.fastq.gz Trimmed/"$1"_R1.trimmed.fastq.gz temp/"$1".matched.index.fq.gz temp/"$1".matched.read1.fq.gz
#/Users/viklj600/BioData/ngsutils/bin/fastqutils properpairs -z -f IndexReads/"$1"_I1.fastq.gz  Trimmed/"$1"_R2.trimmed.fastq.gz temp/"$1".matched.index.fq.gz temp/"$1".matched.read2.fq.gz
