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
from celescope.tools.utils import *
from celescope.tools.Analysis import Analysis
from .otsu import *
import celescope.tools

toolsdir = os.path.dirname(celescope.tools.__file__)

class Analysis_capture_virus(Analysis):

    def get_virus_tsne(self):
        df = pd.merge(self.tsne_df, self.virus_df, on="barcode", how="left")
        self.df_file = f'{self.outdir}/{self.sample}_tsne.tsv'
        df.to_csv(self.df_file, sep='\t')
        df["UMI"] = df["UMI"].fillna(0)
        tSNE_1 = list(df.tSNE_1)
        tSNE_2 = list(df.tSNE_2)
        virus_UMI = list(df.UMI)
        res = {"tSNE_1": tSNE_1, "tSNE_2": tSNE_2, "virus_UMI": virus_UMI}
        return res

    def get_virus_df(self):
        self.virus_df = pd.read_csv(self.virus_file, sep="\t")

    def report(self):
        t = reporter(
            name=self.step,
            assay=self.assay,
            sample=self.sample,
            outdir=self.outdir + '/..',
        )
        t.get_report()

    def run(self):
        cluster_tsne = self.get_cluster_tsne(colname='cluster')
        self.get_virus_df()
        virus_tsne = self.get_virus_tsne()
        table_dict = self.get_marker_gene_table()
        self.report_prepare(
            cluster_tsne=cluster_tsne,
            virus_tsne=virus_tsne,
            table_dict=table_dict,
        )
        self.report()

@add_log
def otsu_thresh(virus_file, outdir, sample):
    df = pd.read_csv(virus_file, sep='\t')
    array = np.log10(df["UMI"])
    hist = array2hist(array)
    thresh = threshold_otsu(hist)
    fname = f'{outdir}/{sample}_otsu_plot.png'
    makePlot(hist, thresh, fname)

    threshold = int(10 ** thresh)
    otsu_thresh.logger.info(f'otsu UMI threshold: {threshold}')
    df_thresh = df[df["UMI"] >= threshold]
    otsu_virus_file = f'{outdir}/{sample}_otsu_count.tsv'
    df_thresh.to_csv(otsu_virus_file, sep='\t')
    return otsu_virus_file


@add_log
def analysis_capture_virus(args):

    # check dir
    outdir = args.outdir
    sample = args.sample
    virus_file = args.virus_file
    match_dir = args.match_dir
    assay = args.assay
    umi_threshold = args.umi_threshold

    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % (outdir))

    if umi_threshold == 'otsu':
        virus_file = otsu_thresh(virus_file, outdir, sample)        

    ana = Analysis_capture_virus(     
        sample,
        outdir,
        assay,
        match_dir=match_dir,
        step='analysis_capture_virus',     
    )
    ana.virus_file = virus_file
    ana.run()


def get_opts_analysis_capture_virus(parser, sub_program):
    parser.add_argument("--umi_threshold", help='method to find virus UMI threshold', 
        choices=['otsu', 'none'], default='otsu')
    if sub_program:
        s_common(parser)
        parser.add_argument('--match_dir', help='match_dir', required=True)
        parser.add_argument('--virus_file', help='virus UMI count file', required=True)
        
