#!/usr/bin/env python3

import io
import os
import pandas as pd
import click

def read_vcf(path):
    lines = []
    with open(path, 'r') as f:
        for l in f:
            # select the CSQ header, assumes "##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format:" - improve.
            if l.startswith('##INFO=<ID=CSQ'):
                csq_holder = l.split(" ")[6].split('|')
            elif not l.startswith('##'):
                lines.append(l)
    return pd.read_table(
        io.StringIO(str.join(os.linesep, lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str}
    ).rename(columns={'#CHROM': 'CHROM'}), csq_holder

def splitdataframe(df, target_column, separator):
    ''' df = dataframe to split,
    target_column = the column containing the values to split
    separator = the symbol used to perform the split
    returns: a dataframe with each entry for the target column separated, with each element moved into a new row.
    The values in the other columns are duplicated across the newly divided rows.
    '''
    def splitListToRows(row, row_accumulator, target_column, separator):
        split_row = row[target_column].split(separator)
        for s in split_row:
            new_row = row.to_dict()
            new_row[target_column] = s
            row_accumulator.append(new_row)
    new_rows = []
    df.apply(splitListToRows, axis=1, args=(
        new_rows, target_column, separator))
    new_df = pd.DataFrame(new_rows)
    return new_df


def expandandextract_techInfo(df, sampleID):
    techinfo = df[sampleID+'.UMI.sorted.trimmed.bam'].str.split(':', expand=True)
    techinfo.columns = ['GT', 'GQ', 'AD', 'DP', 'VF', 'NL', 'SB']
    altDepth = techinfo['AD'].str.split(',', expand=True)
    altDepth.columns = ['REF_DEPTH', 'ALT_DEPTH']
    alldata = pd.concat([df, techinfo, altDepth], axis=1)
    alldataselect = alldata[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'FILTER', 'INFO', 'GQ', 'DP', 'VF', 'REF_DEPTH', 'ALT_DEPTH']]
    return alldataselect

def expandandextract_Info(df, csq, sampleID):
    info = df['INFO'].str.split(';', expand=True)
    info[0] = info[0].str.replace('DP=', '')
    info[1] = info[1].str.replace('AF=', '')
    info[2] = info[2].str.replace('SB=', '')    
    info[3] = info[3].str.replace('DP4=', '')
    info[4] = info[4].str.replace('CSQ=', '')
    info[5] = info[5].str.replace('CNT=', '')
    info.columns = ['DP', 'AF', 'SB', 'DP4', 'CSQ', 'CNT']
    alldata = pd.concat([df, info], axis=1)
    expandedcsq = splitdataframe(alldata, 'CSQ', ',')
    csqextracted = expandedcsq["CSQ"].str.split('|', expand=True)
    csqextracted.columns = csq
    alldata_csq = pd.concat([expandedcsq, csqextracted], axis=1)
    alldata_csq['SAMPLE']=sampleID
    alldata_csq_select = alldata_csq[['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'SYMBOL', 'Consequence', 'BIOTYPE', 'VARIANT_CLASS', 'Feature_type', 'IMPACT', 'FLAGS', 'DP', 'AF', 'Existing_variation', 'EUR_AF', 'gnomAD_NFE_AF', 'MAX_AF', 'ID', 'CNT', 'FILTER', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons']]
    alldata_csq_renamed = alldata_csq_select.rename(columns={'Consequence': 'CONSEQUENCE', 'Feature_type': 'FEATURE_TYPE', 'DP': 'DEPTH', 'VF': 'VAR', 'ID': 'COSMIC_IDs', 'CNT': 'COSMIC_CASES'})
    return alldata_csq_renamed



@click.command()
@click.option(
    '--inputfile', '-i',
    help='VEP VCF with added Cosmic info',
)
@click.option(
    '--outputfile', '-o',
    help='outputfile!',
)
@click.option(
    '--samplename', '-s',
    help='Sample name',
)

def main(inputfile, outputfile, samplename):
    testdb, csq = read_vcf(inputfile)
    #testdb_split_techinfo = expandandextract_techInfo(testdb, samplename)
    testdb_split_info = expandandextract_Info(testdb, csq, samplename)
    testdb_split_info.to_csv('{}'.format(outputfile), sep='\t', index=False)

if __name__ == '__main__':
    main()
