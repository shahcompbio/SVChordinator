import os
import re
import sys
import pandas as pd
import argparse
import numpy as np


def parse_args():
    desc = "convert minda VCF to tsv and "
    p = argparse.ArgumentParser(description=desc)
    p.add_argument('-input', '--input', help="input sniffles vcf", required=True)
    p.add_argument('-output', '--output', help="output tsv", required=True)
    args = p.parse_args()
    return args


def proc_info(row):
    field = row['INFO'].split(';')
    infos = {}
    for each_info in field:
        if '=' in each_info:
            key, value = each_info.split('=')
            infos[key] = value  # infos['SVTYPE'] -> 'BND'
        else:
            infos[each_info] = True  # infos['IMPRECISE'] -> True
    return infos


def get_chrom2_pos2(row, infos):
    assert 'SVTYPE' in infos, f'infos does not have SVTYPE:\n{row}'
    # get for translocation
    if infos['SVTYPE'] == 'BND' or infos['SVTYPE'] == 'TRA':  # ALT: N]chrX:2212928]
        pat = re.search('[\[\]](.+):(\d+)', row['ALT'])
        assert pat
        chrom2, pos2 = pat.groups()
        print(chrom2)
    elif infos['SVTYPE'] == 'INS':
        chrom2, pos2 = row['#CHROM'], row['POS']
    else:
        chrom2 = row['#CHROM']
        pos2 = int(row['POS']) + int(infos['SVLEN'])
    return chrom2, pos2


def get_svlen(infos):
        return int(infos['SVLEN'])


class MAF:
    def __init__(self, tumor_id):
        self.chrom1s = []
        self.pos1s = []
        self.chrom2s = []
        self.pos2s = []
        self.refs = []
        self.alts = []
        self.filters = []
        self.svtypes = []
        self.svlens = []
        self.callers = []
        self.ID = []
        self.strand1s = []
        self.strand2s = []
        self.tumor_id = tumor_id

    def __repr__(self):
        return f'MAF of tumor_id {self.tumor_id}'

    def proc_vcf(self, vcf_path):
        if vcf_path.endswith('.gz'):
            vcf_file = gzip.open(vcf_path, 'rt', encoding='utf-8')
        else:
            vcf_file = open(vcf_path, 'r')
        for line in vcf_file:
            if line.startswith('##'):
                continue
            elif line.startswith('#'):
                header = line.rstrip().split('\t')
                assert len(header) <= 11, f'ERROR: header\n{header}\ntoo long'
                for ix in range(9, len(header)):
                    if header[ix] == self.tumor_id:  # tumor
                        header[ix] = 'TUMOR'
                    ## we can fix this when we move to somatic calls
                    # else:
                    #     header[ix] = 'NORMAL'
                self.header = header
                continue
            field = line.strip().split('\t')
            row = dict(zip(self.header, field))

            self.add_data(row)

    def add_data(self, row):
        infos = proc_info(row)
        remove_sv = False
        chrom2, pos2 = get_chrom2_pos2(row, infos)
        svlen = get_svlen(infos)
        assert len(self.info_keys & set(infos.keys())) == len(self.info_keys)
        if len(infos['STRANDS']) == 2:
            strand1, strand2 = infos['STRANDS']
        else:
            strand1 = strand2 = infos['STRANDS']

        if not remove_sv:
            self.chrom1s.append(row['#CHROM'])
            self.pos1s.append(int(row['POS']))
            self.chrom2s.append(chrom2)  # BND?
            self.pos2s.append(pos2)
            self.ID.append(row['ID'])
            self.refs.append(row['REF'])
            self.alts.append(row['ALT'])
            self.filters.append(row['FILTER'])
            self.svtypes.append(infos['SVTYPE'])
            self.svlens.append(svlen)
            self.callers.append(infos["SUPP_VEC"])
            self.strand1s.append(strand1)
            self.strand2s.append(strand2)
## for now, I'm giving this the same format as the fusions and nanomonsv
    def to_df(self):
        self.df = pd.DataFrame({
            'ID': self.ID,
            'chrom1': self.chrom1s,
            'base1': self.pos1s,
            'strand1': self.strand1s,
            'chrom2': self.chrom2s,
            'base2': self.pos2s,
            'strand2': self.strand2s,
            'Ref_Seq': self.refs,
            'Alt_Seq': self.alts,
            'SV_Type': self.svtypes,
            'SV_Len': self.svlens,
            'SV_callers': self.callers,
        })

     #   self.df['rnames'] = self.rnames

        return self.df


if __name__ == "__main__":
    # args = parse_args()
    mafobj = MAF('SAMPLE')
    vcf_path = '/Users/asherpreskasteinberg/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/lab_notebook/APS017.1_3x3_SV_analysis/variants2table/SHAH_H002190_T02_01_WG01-T.nanomonsv.result.vcf'
    mafobj.proc_vcf(vcf_path)
    sniffles_maf = mafobj.to_df()
    maf = mafobj.to_df()
    #maf.to_csv(snakemake.output["tsv"], sep='\t', index=False)
    csv_path = '/Users/asherpreskasteinberg/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/lab_notebook/APS017.1_3x3_SV_analysis/variants2table/SHAH_H002190_T02_01_WG01-T.nanomonsv.result.tsv'
    maf.to_csv(csv_path, sep='\t', index=False)
