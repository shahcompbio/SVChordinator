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
    p.add_argument("-read_support", "--read_support", type=bool, default=False,
                   help="pull read support info from sniffles vcf; \
                   if True, --sniffles is required")
    p.add_argument("-sniffles_vcf", "--sniffles_vcf",
                   help="sniffles vcf for calculating read support; required if \
                   --read_support flag is True", required=False)
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


def proc_samples(row):
    formats = row['FORMAT'].split(':')
    tumor_gts = row['TUMOR'].split(':')
    assert 'DR' in formats
    has_normal = False
    if 'NORMAL' in row:
        has_normal = True
        normal_gts = row['NORMAL'].split(':')
        return (has_normal,
                (dict(zip(formats, tumor_gts)),
                 dict(zip(formats, normal_gts))))
    return (has_normal,
            dict(zip(formats, tumor_gts)))

## for translocations
def get_chrom2_pos2(row, infos):
    assert 'SVTYPE' in infos, f'infos does not have SVTYPE:\n{row}'
    if infos['SVTYPE'] == 'BND' or infos['SVTYPE'] == 'TRA':  # ALT: N]chrX:2212928]
        pat = re.search('[\[\]](.+):(\d+)', row['ALT'])
        assert pat
        chrom2, pos2 = pat.groups()
        print(chrom2)
    elif 'END' in infos:
        chrom2, pos2 = row['#CHROM'], infos['END']
    else:
        print(f'row:\n{row}\ninfos:\n{infos}')
    return chrom2, int(pos2)


def get_svlen(infos):
    if infos['SVTYPE'] == 'BND':
        return 3e9
    else:
        return abs(int(infos['SVLEN']))


class MAF:
    def __init__(self, tumor_id, survivor=True):
        self.chrom1s = []
        self.pos1s = []
        self.strand1s = []
        self.chrom2s = []
        self.pos2s = []
        self.strand2s = []
        self.refs = []
        self.alts = []
        self.filters = []
        self.svtypes = []
        self.svlens = []
        self.tumor_DRs = []
        self.tumor_DVs = []
        self.tumor_VAFs = []
        self.normal_DRs = []
        self.normal_DVs = []
        self.normal_VAFs = []
        self.rnames = []
        self.sniffles_ID = []
        if survivor:
            self.info_keys = {'SVTYPE', 'STRANDS', }
        else:
            self.info_keys = {'SVTYPE', 'STRAND', }
        self.tumor_id = tumor_id  # tumor id only
        self.has_normal = False  # does GT include normal

    def __repr__(self):
        return f'MAF of tumor_id {self.tumor_id}'

    def proc_vcf(self, vcf_path, survivor=True):
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

            self.add_data(row, survivor=survivor)

    def add_data(self, row, survivor=True):
        infos = proc_info(row)
        self.has_normal, samples = proc_samples(row)
        if self.has_normal:
            tumors, normals = samples
        else:
            tumors = samples
        remove_sv = False
        # if tumors['GT'] in ('./.', '.|.', '0/0', '0|0',):
        #    remove_sv = True
        chrom2, pos2 = get_chrom2_pos2(row, infos)
        svlen = get_svlen(infos)
        # rnames = infos['RNAMES']  # $rname1,$rname2,...

        assert len(self.info_keys & set(infos.keys())) == len(self.info_keys)
        if survivor:
            if len(infos['STRANDS']) == 2:
                strand1, strand2 = infos['STRANDS']
            else:
                strand1 = strand2 = infos['STRANDS']
        else:
            if len(infos['STRAND']) == 2:
                strand1, strand2 = infos['STRAND']
            else:
                strand1 = strand2 = infos['STRAND']
        if not remove_sv:
            self.chrom1s.append(row['#CHROM'])
            self.pos1s.append(int(row['POS']))
            self.strand1s.append(strand1)
            self.chrom2s.append(chrom2)  # BND?
            self.pos2s.append(pos2)
            self.strand2s.append(strand2)
            self.refs.append(row['REF'])
            self.alts.append(row['ALT'])
            self.filters.append(row['FILTER'])
            self.svtypes.append(infos['SVTYPE'])
            self.svlens.append(svlen)
            # self.rnames.append(rnames)
            ## split the DR value as this has a funny format in survivor
            if survivor:
                self.tumor_DRs.append(np.nan)
                self.tumor_DVs.append(np.nan)
                self.tumor_VAFs.append(np.nan)
                ## also add sniffles ID ....
                sniffles_line = row['TUMOR']
                terms = sniffles_line.split(":")
                self.sniffles_ID.append(terms[7])
            else:
                _, _, DR, DV = row["TUMOR"].split(":")
                self.tumor_DRs.append(int(DR))
                self.tumor_DVs.append(int(DV))
                #self.tumor_VAFs.append(infos["AF"])
                self.sniffles_ID.append(row["ID"])
                ## compute AF
                if (float(DV) + float(DR)) > 0:
                    tumor_VAF = float(DV) / (float(DV) + float(DR))
                else:
                    tumor_VAF = 0.0
                self.tumor_VAFs.append(tumor_VAF)

            if self.has_normal:
                self.normal_DRs.append(int(normals['DR']))
                self.normal_DVs.append(int(normals['DV']))
                if (float(normals['DV']) + float(normals['DR'])) > 0:
                    normal_VAF = float(normals['DV']) / (float(normals['DV']) + float(normals['DR']))
                else:
                    normal_VAF = 0.0
                self.normal_VAFs.append(normal_VAF)
## for now, I'm giving this the same format as the fusions and nanomonsv
    def to_df(self):
        self.df = pd.DataFrame({
            'chrom1': self.chrom1s,
            'base1': self.pos1s,
            'sniffles_ID': self.sniffles_ID,
            'strand1': self.strand1s,
            'chrom2': self.chrom2s,
            'base2': self.pos2s,
            'strand2': self.strand2s,
            'Ref_Seq': self.refs,
            'Alt_Seq': self.alts,
            'FILTER': self.filters,
            'SV_Type': self.svtypes,
            'SV_LEN': self.svlens,
            'Supporting_Read_Num_Ref': self.tumor_DRs,  # reference reads
            'Supporting_Read_Num_Tumor': self.tumor_DVs,  # variant reads
            'tumor_VAF': self.tumor_VAFs,
            # actually means VAF: https://github.com/fritzsedlazeck/Sniffles/blob/313c6001421738afafcec6dbf4dd24e41223a668/src/sniffles/postprocessing.py#L273
        })
        if self.has_normal:
            self.df['normal_DR'] = self.normal_DRs  # reference reads
            self.df['normal_DV'] = self.normal_DVs  # variant reads
            self.df['normal_VAF'] = self.normal_VAFs  # VAFs

     #   self.df['rnames'] = self.rnames

        return self.df


if __name__ == "__main__":
    # args = parse_args()
    mafobj = MAF('SAMPLE', survivor=False)
    # vcf_path = args.sniffles_vcf
    vcf_path = snakemake.input["sniffles_vcf"]
    mafobj.proc_vcf(vcf_path, survivor=False)
    sniffles_maf = mafobj.to_df()
    maf = mafobj.to_df()
    maf.to_csv(snakemake.output["tsv"], sep='\t', index=False)
    ## split into many mafs for gene annotation step
    os.makedirs(snakemake.params["split_out"], exist_ok=True)
    out_paths = snakemake.output["split_tsv"]
    num_splits = len(out_paths)
    # Split the DataFrame
    split_dfs = np.array_split(maf, num_splits)
    # Convert split arrays back into DataFrames
    split_dfs = [pd.DataFrame(split) for split in split_dfs]
    for i in np.arange(0, len(split_dfs)):
        split_df = split_dfs[i]
        split_path = out_paths[i]
        split_df.to_csv(split_path, sep='\t', index=False)
    # with open(args.rnames, 'w') as out:
    #     for rnames in maf['rnames'].values:
    #         lines = rnames.replace(',', '\n') + '\n'
    #         out.write(lines)
