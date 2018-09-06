process coverage_analysis_standard{
    tag {idSample}

    publishDir directoryMap.coverage, mode: 'link'

    input:
      set idPatient, idSample, file(bam), file(bai) from trimmed_StandardBAM
      file(regions) from Channel.value(regionsFile)

    output:
      file("${idSample}.standard.coverage.txt") into standard_coverage

    script:
    """
    sambamba depth region --filter 'mapping_quality > 20' -q 20 -T 10 -T 50 -T 100 -T 500 -T 1000 -L ${regions} -o ${idSample}.rawcoverage.txt ${bam}
    ParseSamCoverage.py -i ${idSample}.rawcoverage.txt -o ${idSample}.standard.coverage.txt
    """
}

process coverage_analysis_UMI{
    tag {idSample}

    publishDir directoryMap.coverage, mode: 'link'

    input:
      set idPatient, idSample, file(bam), file(bai) from trimmed_umiBAM_forcoverage
      file(regions) from Channel.value(regionsFile)

    output:
      file("${idSample}.UMI.coverage.txt") into umi_coverage

    script:
    """
    sambamba depth region --filter 'mapping_quality > 20' -q 20 -T 10 -T 50 -T 100 -T 500 -T 1000 -L ${regions} -o ${idSample}.rawcoverage.txt ${bam}
    ParseSamCoverage.py -i ${idSample}.rawcoverage.txt -o ${idSample}.UMI.coverage.txt
    """
}