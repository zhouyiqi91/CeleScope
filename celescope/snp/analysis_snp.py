import os
import sys
import json
import logging
import re
import numpy as np
import pandas as pd
import glob
from scipy.io import mmwrite
from scipy.sparse import csr_matrix
from celescope.tools.report import reporter
from celescope.tools.utils import glob_genomeDir, log
from celescope.tools.utils import cluster_tsne_list, marker_table, report_prepare, parse_vcf
from celescope.snp.snpCalling import read_index
import celescope.tools

toolsdir = os.path.dirname(celescope.tools.__file__)

def vcf_data(vcf_file, tsne_df, index_file):
    """
    return data dic
    """
    # df_vcf_count
    df_vcf = parse_vcf(vcf_file)
    df_vcf.index = [int(index) for index in list(df_vcf['CELL'])]
    df_vcf_count = df_vcf.groupby(['CHROM','POS','ALLELES','GENE']).agg({
        'CELL':'count'}).reset_index().sort_values(['CELL'],ascending=False)
    df_vcf_count_top10 = df_vcf_count.iloc[0:9]

    # vcf_tsne
    df_index, df_valid = read_index(index_file)
    df_vcf_barcode = pd.merge(df_vcf, df_valid, left_index=True, right_index=True, how="left")
    df_vcf_tsne = pd.merge(df_vcf_barcode, tsne_df, on="barcode", how="left")

    # vcf_tsne_all
    df_vcf_tsne_count = df_vcf_tsne.groupby(['barcode','CELL']).agg({'ALLELES':'count'})
    df_vcf_tsne_all = pd.merge(tsne_df, df_vcf_tsne_count, on='barcode', how='left')
    df_index = df_index.reset_index()
    df_vcf_tsne_all = pd.merge(df_vcf_tsne_all, df_index, on='barcode', how='left')
    df_vcf_tsne_all.fillna(0, inplace=True)
    print(df_vcf_tsne_all)
    tSNE_1 = list(df_vcf_tsne_all.tSNE_1)
    tSNE_2 = list(df_vcf_tsne_all.tSNE_2)
    def return_text(row):
        text = f'Variants:{str(int(row["ALLELES"]))}'
        return text
    text = list(df_vcf_tsne_all.apply(return_text, axis=1))
    value = list(df_vcf_tsne_all.ALLELES)
    print(text)
    print(value)
    title = 't-SNE plot Colored by Cell Variant Counts'
    res = {"tSNE_1": tSNE_1, "tSNE_2": tSNE_2, "text": text, 'value':value, 'title':title}
    return res


@log
def analysis_snp(args):

    outdir = args.outdir
    sample = args.sample
    match_dir = args.match_dir
    vcf_file = args.vcf_anno
    index_file = args.index_file

    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % (outdir))

    # report
    tsne_df_file = glob.glob(f'{match_dir}/*analysis*/*tsne_coord.tsv')[0]
    marker_df_file = glob.glob(f'{match_dir}/*analysis*/*markers.tsv')[0]
    tsne_df = pd.read_csv(tsne_df_file, sep="\t")
    tsne_df.rename(columns={"Unnamed: 0": "barcode"}, inplace=True)
    marker_df = pd.read_csv(marker_df_file, sep="\t")
    #virus_df = pd.read_csv(virus_file, sep="\t")

    cluster_tsne = cluster_tsne_list(tsne_df)
    count_tsne = vcf_data(vcf_file, tsne_df, index_file)

    report_prepare(outdir, cluster_tsne=cluster_tsne, count_tsne=count_tsne)

    t = reporter(
        name='analysis_snp',
        assay=args.assay,
        sample=args.sample,
        outdir=args.outdir + '/..')
    t.get_report()


def get_opts_analysis_snp(parser, sub_program):
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument('--match_dir', help='match_dir', required=True)
        parser.add_argument('--vcf_anno', help='annotated vcf file', required=True)
        parser.add_argument('--index_file', help='index_file', required=True)
        parser.add_argument('--assay', help='assay', required=True)