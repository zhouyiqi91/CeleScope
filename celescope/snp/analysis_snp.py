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

class analysis_variant():

    def __init__(self, outdir, sample, match_dir, vcf_file, index_file, assay):
        self.outdir = outdir
        self.sample = sample
        self.match_dir = match_dir
        self.vcf_file = vcf_file
        self.index_file = index_file 
        self.assay = assay  
        if not os.path.exists(outdir):
            os.system('mkdir -p %s' % (outdir))

        #
        tsne_df_file = glob.glob(f'{match_dir}/*analysis*/*tsne_coord.tsv')[0]
        #marker_df_file = glob.glob(f'{match_dir}/*analysis*/*markers.tsv')[0]
        self.tsne_df = pd.read_csv(tsne_df_file, sep="\t")
        self.tsne_df.rename(columns={"Unnamed: 0": "barcode"}, inplace=True)

        self.cluster_tsne = cluster_tsne_list(self.tsne_df)

    def vcf_data(self):
        """
        return data dic
        """
        # df_vcf_count
        df_vcf = parse_vcf(self.vcf_file)
        df_vcf.index = [int(index) for index in list(df_vcf['CELL'])]
        df_vcf_count = df_vcf.groupby(['CHROM','POS','ALLELES','GENE']).agg({
            'CELL':'count'}).reset_index().sort_values(['CELL'],ascending=False)

        # vcf_tsne
        df_index, df_valid = read_index(self.index_file)
        df_vcf_barcode = pd.merge(df_vcf, df_valid, left_index=True, right_index=True, how="left")
        df_vcf_tsne = pd.merge(df_vcf_barcode, self.tsne_df, on="barcode", how="left")

        # df_table
        df_vcf_count.rename(columns={'CELL':'NCELL'}, inplace=True)
        df_vcf_tsne.rename(columns={'CELL':'INDEX'}, inplace=True)
        df_vcf_tsne = df_vcf_tsne.astype({'CHROM':'object', 'POS':'int64'})
        df_vcf_count = df_vcf_count.astype({'CHROM':'object', 'POS':'int64'})
        df_table = pd.merge(df_vcf_tsne, df_vcf_count, on=['CHROM',"POS","ALLELES",'GENE'], how='left')
        self.df_table_full = df_table
        df_table_full_file = f'{self.outdir}/{self.sample}_variant_table.tsv'
        self.df_table_full.to_csv(df_table_full_file, sep='\t')
        cols = ['CHROM', 'POS', 'GENE', 'ALLELES', 'DP', 'GT', 'INDEX', 'cluster', 'NCELL']
        df_table = df_table.loc[:, cols]
        df_table.rename(columns={'cluster':'CLUSTER'}, inplace=True)
        self.df_table = df_table

        # vcf_tsne_all
        df_vcf_tsne_count = df_vcf_tsne.groupby(['barcode','INDEX']).agg({'ALLELES':'count'})
        df_vcf_tsne_all = pd.merge(self.tsne_df, df_vcf_tsne_count, on='barcode', how='left')
        df_index = df_index.reset_index()
        df_vcf_tsne_all = pd.merge(df_vcf_tsne_all, df_index, on='barcode', how='left')
        df_vcf_tsne_all.fillna(0, inplace=True)
        tSNE_1 = list(df_vcf_tsne_all.tSNE_1)
        tSNE_2 = list(df_vcf_tsne_all.tSNE_2)
        def return_text(row):
            text = f'Variants:{str(int(row["ALLELES"]))}<br>Cell_Index:{row["cell_index"]}'
            return text
        text = list(df_vcf_tsne_all.apply(return_text, axis=1))
        value = list(df_vcf_tsne_all.ALLELES)
        title = 't-SNE plot Colored by Cell Variant Counts'
        self.count_tsne = {"tSNE_1": tSNE_1, "tSNE_2": tSNE_2, "text": text, 'value':value, 'title':title}

    def get_variant_table(self):
        """
        return html code
        """
        self.variant_table = {}
        self.variant_table['title'] = 'Variant Table'
        self.variant_table['table'] = self.df_table.to_html(
            escape=False,
            index=False,
            table_id="variant_table",
            justify="center")

    def report(self):
        t = reporter(
        name='analysis_snp',
        assay=self.assay,
        sample=self.sample,
        outdir=self.outdir + '/..')
        t.get_report()

    def run(self):
        self.vcf_data()
        self.get_variant_table()
        report_prepare(
            self.outdir,
            cluster_tsne=self.cluster_tsne, 
            count_tsne=self.count_tsne, 
            variant_table=self.variant_table
        )
        self.report()

@log
def analysis_snp(args):
    step_snp = analysis_variant(
        args.outdir,
        args.sample,
        args.match_dir,
        args.vcf_anno,
        args.index_file,
        args.assay
    )
    step_snp.run()

def get_opts_analysis_snp(parser, sub_program):
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument('--match_dir', help='match_dir', required=True)
        parser.add_argument('--vcf_anno', help='annotated vcf file', required=True)
        parser.add_argument('--index_file', help='index_file', required=True)
        parser.add_argument('--assay', help='assay', required=True)