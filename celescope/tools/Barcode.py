#!/bin/env python
# coding=utf8
import os
import re
import io
import gzip
import subprocess
import sys
import glob
import pandas as pd
import pysam
from collections import defaultdict, Counter
from itertools import combinations, permutations, islice
from xopen import xopen
from celescope.tools.utils import format_number, log, seq_ranges, read_fasta, genDict, read_one_col
from celescope.tools.report import reporter
from celescope.tools.__init__ import __PATTERN_DICT__

class Barcode():
    def __init__(self, fq1, fq2, sample, outdir): 
        # TODO
        self.read = 0
        self.nRead = 10000
        self.T4 = 0
        self.L57C = 0
        self.fq1_file, self.fq2_file = self.fq1, self.fq2

    @log
    def merge_fastq(self):
        '''
        merge fastq if len(fq1) > 1
        '''
        fq1_list = self.fq1.split(",")
        fq2_list = self.fq2.split(",")
        if len(fq1_list) == 0:
            raise Exception('empty fastq file path!')
        elif len(fq1_list) == 1:
            self.fq1_file, self.fq2_file = self.fq1, self.fq2
        else:
            fastq_dir = f'{self.outdir}/../merge_fastq'
            if not os.path.exists(fastq_dir):
                os.system('mkdir -p %s' % fastq_dir)
            fq1_file = f"{fastq_dir}/{self.sample}_1.fq.gz"
            fq2_file = f"{fastq_dir}/{self.sample}_2.fq.gz"
            fq1_files = " ".join(fq1_list)
            fq2_files = " ".join(fq2_list)
            fq1_cmd = f"cat {fq1_files} > {fq1_file}"
            fq2_cmd = f"cat {fq2_files} > {fq2_file}"
            Barcode.merge_fastq.logger.info(fq1_cmd)
            os.system(fq1_cmd)
            Barcode.merge_fastq.logger.info(fq2_cmd)
            os.system(fq2_cmd)
            self.fq1_file = fq1_file
            self.fq2_file = fq2_file
        return self.fq1_file, self.fq2_file

    '''
    @log
    def get_chemistry(self):
        '''
        'scopeV2.0.1': 'C8L16C8L16C8L1U8T18'
        'scopeV2.1.1': 'C8L16C8L16C8L1U12T18'
        'scopeV2.2.1': 'C8L16C8L16C8L1U12T18' with 4 types of linkers
        '''
        # init
        linker_4_file, _whitelist = Barcode.get_scope_bc('scopeV2.2.1')
        linker_4_list, _num = read_one_col(linker_4_file)
        linker_4_dict = defaultdict(int)
        pattern_dict = parse_pattern('C8L16C8L16C8L1U12T18')

        with pysam.FastxFile(self.fq1) as fh:
            for i in range(self.nRead):
                entry = fh.__next__()
                seq = entry.sequence
                L57C = seq[56]
                if L57C == 'C':
                    self.L57C += 1
                T4 = seq[65:69]
                if T4 == 'TTTT':
                    self.T4 += 1
                linker = seq_ranges(seq, pattern_dict=pattern_dict)
                if linker in linker_4_list:
                    linker_4_dict[linker] += 1

        self.percent_T4 =  self.T4 / self.nRead 
        self.percent_L57C = self.L57C / self.nRead
        Chemistry.get_chemistry.logger.info(f'percent T4: {self.percent_T4}')
        Chemistry.get_chemistry.logger.info(f'percent L57C: {self.percent_L57C}')
        if self.percent_T4 > 0.5:
            self.chemistry = 'scopeV2.0.1'
        else:
            # V2.1.1 or V2.2.1 or failed
            self.valid_linker_type = 0
            for linker in linker_4_dict:
                linker_4_dict[linker] = linker_4_dict[linker] / self.nRead
                if linker_4_dict[linker] > 0.05:
                    self.valid_linker_type += 1
            Chemistry.get_chemistry.logger.info(self.linker_4_dict)
            if self.valid_linker_type == 0:
                raise Exception('auto chemistry detection failed!')
            elif self.valid_linker_type == 1:
                self.chemistry == 'scopeV2.1.1'
            elif self.valid_linker_type < 4:
                self.chemistry == 'scopeV2.1.1'
                Chemistry.get_chemistry.logger.warning(
                    f'chemistry scopeV2.2.1 only has {self.valid_linker_type} linker types!')
            else:
                self.chemistry == 'scopeV2.2.1'
        Chemistry.get_chemistry.logger.info(f'chemistry: {self.chemistry}')
        return self.chemistry
    '''

    @staticmethod
    def get_scope_bc(bctype):
        import celescope
        root_path = os.path.dirname(celescope.__file__)
        linker_f = glob.glob(f'{root_path}/data/chemistry/{bctype}/linker*')[0]
        whitelist_f = f'{root_path}/data/chemistry/{bctype}/bclist'
        return linker_f, whitelist_f