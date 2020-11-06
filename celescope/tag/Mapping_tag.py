from celescope.tools.report import reporter
from celescope.tools.utils import format_number, log, genDict
from celescope.tools.utils import read_fasta, process_read, gen_stat
from celescope.tools.barcode import parse_pattern
import gzip
import os
import pandas as pd
import logging
import numpy as np
import sys
import argparse
import re
import glob
from collections import defaultdict
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


class Mapping_tag():

    def __init__(
        self,
        sample,
        outdir,
        assay,
        fq,
        fq_pattern,
        linker_fasta,
        barcode_fasta,
        ):

        # read args
        self.sample = sample
        self.outdir = outdir
        self.assay = assay
        self.fq = fq
        self.fq_pattern = fq_pattern
        self.linker_fasta = linker_fasta
        self.barcode_fasta = barcode_fasta

        # process
        self.barcode_dic, self.barcode_length = read_fasta(self.barcode_fasta, equal=True)
        if self.linker_fasta and self.linker_fasta != 'None':
            self.linker_dic, self.linker_length = read_fasta(self.linker_fasta, equal=True)
        else:
            self.linker_dic, self.linker_length = {}, 0
        self.pattern_dict = parse_pattern(self.fq_pattern)

        # check barcode length
        barcode1 = self.pattern_dict["C"][0]
        # end - start
        pattern_barcode_length = barcode1[1] - barcode1[0]
        if pattern_barcode_length != self.barcode_length:
            raise Exception(
                f'''barcode fasta length {self.barcode_length} 
                != pattern barcode length {pattern_barcode_length}'''
            )

        self.res_dic = genDict()
        self.res_sum_dic = genDict(dim=2)

        self.match_barcode = []
        self.read_count_file = f'{outdir}/{sample}_read_count.tsv'
        self.UMI_count_file = f'{outdir}/{sample}_UMI_count.tsv'
        self.stat_file = self.outdir + "/stat.txt"

        if not os.path.exists(outdir):
            os.system('mkdir -p %s' % outdir)

    def read_to_dic(self):

        self.res_dic, self.metrics = process_read(
            self.fq,
            self.pattern_dict,
            self.barcode_dic,
            self.linker_dic,
            self.barcode_length,
            self.linker_length,
        )

        # write dic to pandas df
        rows = []
        for barcode in self.res_dic:
            for tag_name in self.res_dic[barcode]:
                for umi in self.res_dic[barcode][tag_name]:
                    rows.append([barcode, tag_name, umi,
                                 self.res_dic[barcode][tag_name][umi]])
        df_read_count = pd.DataFrame(rows)
        df_read_count.rename(
            columns={
                0: "barcode",
                1: "tag_name",
                2: "UMI",
                3: "read_count"},
            inplace=True)
        df_read_count.to_csv(
            self.read_count_file, sep="\t", index=False)

        # stat
        stat = pd.DataFrame(
            {
                "item": [
                    "Reads Mapped",
                    "Reads Unmapped too Short",
                    "Reads Unmapped Invalid Linker",
                    "Reads Unmapped Invalid Barcode",
                ],
                "count": [
                    self.metrics["Reads Mapped"],
                    self.metrics["Reads Unmapped too Short"],
                    self.metrics["Reads Unmapped Invalid Linker"],
                    self.metrics["Reads Unmapped Invalid Barcode"],
                ],
            }, 
            columns=["item", "count"]
        )
        stat['total_count'] = self.metrics["Total Reads"]
        gen_stat(stat, self.stat_file)

    def tag_count(self):
        for barcode in self.res_dic:
            for tag_name in self.res_dic[barcode]:
                self.res_sum_dic[barcode][tag_name] = len(
                    self.res_dic[barcode][tag_name])

        df_umi_count = pd.DataFrame(self.res_sum_dic)
        df_umi_count = df_umi_count.T
        df_umi_count.fillna(0, inplace=True)
        df_umi_count = df_umi_count.astype(int)
        df_umi_count.to_csv(self.UMI_count_file, sep="\t")

    def report(self):
        t = reporter(
            name='mapping_tag', 
            assay=self.assay, 
            sample=self.sample,
            stat_file=self.stat_file, 
            outdir=self.outdir + '/..')
        t.get_report()

    @log
    def run(self):
        self.read_to_dic()
        self.tag_count()
        self.report()



