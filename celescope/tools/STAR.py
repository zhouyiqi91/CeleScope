#!/bin/env python
# coding=utf8

import os
import re
import sys
import json
import logging
import subprocess
import glob
import pandas as pd
from celescope.tools.utils import format_number, log, format_stat
from celescope.tools.utils import glob_genomeDir
from celescope.tools.report import reporter

class Step_mapping():

    def __init__(self, sample, outdir, assay, thread, fq, genomeDir, out_unmapped, debug, outFilterMatchNmin):
        self.sample = sample
        self.outdir = outdir
        self.assay = assay
        self.thread = thread
        self.fq = fq
        self.genomeDir = genomeDir
        self.out_unmapped = out_unmapped
        self.debug = debug
        self.outFilterMatchNmin = outFilterMatchNmin
        self.refFlat, self.gtf = glob_genomeDir(genomeDir)
        if not os.path.exists(outdir):
            os.system('mkdir -p %s' % (outdir))
        self.stats = pd.Series()
        self.stats_file = f'{outdir}/stat.txt'

    def format_stat(self):
        fh1 = open(self.STAR_map_log, 'r')
        UNIQUE_READS = []
        MULTI_MAPPING_READS = []
        for line in fh1:
            if line.strip() == '':
                continue
            if re.search(r'Uniquely mapped reads', line):
                UNIQUE_READS.append(line.strip().split()[-1])
            if re.search(r'of reads mapped to too many loci', line):
                MULTI_MAPPING_READS.append(line.strip().split()[-1])
        fh1.close()

        fh2 = open(self.picard_region_log, 'r')
        region_dict = {}
        while True:
            line = fh2.readline()
            if not line:
                break
            if line.startswith('## METRICS CLASS'):
                header = fh2.readline().strip().split('\t')
                data = fh2.readline().strip().split('\t')
                region_dict = dict(zip(header, data))
                break
        fh2.close()
        
        Total = float(region_dict['PF_ALIGNED_BASES'])
        Exonic_Regions = int(region_dict['UTR_BASES']) + \
            int(region_dict['CODING_BASES'])
        Intronic_Regions = int(region_dict['INTRONIC_BASES'])
        Intergenic_Regions = int(region_dict['INTERGENIC_BASES'])

        region_dict['Exonic_Regions'] = "{}({:.2%})".format(
            format_number(Exonic_Regions), Exonic_Regions / Total)
        region_dict['Intronic_Regions'] = "{}({:.2%})".format(
            format_number(Intronic_Regions), Intronic_Regions / Total)
        region_dict['Intergenic_Regions'] = "{}({:.2%})".format(
            format_number(Intergenic_Regions), Intergenic_Regions / Total)

        self.stats = self.stats.append(pd.Series(
            f'{format_number(int(UNIQUE_READS[0]))}({UNIQUE_READS[1]})',
            index=['Uniquely Mapped Reads']
        ))
        self.stats = self.stats.append(pd.Series(
            f'{format_number(int(MULTI_MAPPING_READS[0]))}({MULTI_MAPPING_READS[1]})',
            index=['Multi-Mapped Reads']
        ))
        # ribo
        if self.debug:
            f = open(self.ribo_log, 'r')
            for line in f:
                if line.find('#Matched') != -1:
                    items = line.split()
                    Reads_Mapped_to_rRNA = int(items[1])
                if line.find('#Total') != -1:
                    items = line.split()
                    Reads_Total = int(items[1])

            self.stats = self.stats.append(pd.Series(
                format_stat(Reads_Mapped_to_rRNA, Reads_Total),
                index=['Reads Mapped to rRNA']
            ))
            f.close()

        self.stats = self.stats.append(pd.Series(
            region_dict['Exonic_Regions'],
            index=['Base Pairs Mapped to Exonic Regions']
        ))
        self.stats = self.stats.append(pd.Series(
            region_dict['Intronic_Regions'],
            index=['Base Pairs Mapped to Intronic Regions']
        ))
        self.stats = self.stats.append(pd.Series(
            region_dict['Intergenic_Regions'],
            index=['Base Pairs Mapped to Intergenic Regions']
        ))
        self.plot = {'region_labels': ['Exonic Regions', 'Intronic Regions', 'Intergenic Regions'],
                'region_values': [Exonic_Regions, Intronic_Regions, Intergenic_Regions]}   


        self.stats.to_csv(self.stats_file, sep=':', header=False)

    @log
    def ribo(self):
        import celescope
        root_path = os.path.dirname(celescope.__file__)
        human_ribo_fa = f'{root_path}/data/rRNA/human_ribo.fasta'
        self.ribo_log = f'{self.outdir}/{self.sample}_ribo_log.txt'
        self.ribo_run_log = f'{self.outdir}/{self.sample}_ribo_run.log'
        cmd = (
            f'bbduk.sh '
            f'in1={self.fq} '
            f'ref={human_ribo_fa} '
            f'stats={self.ribo_log} '
            f'overwrite=t '
            f'> {self.ribo_run_log} 2>&1 '
        )
        Step_mapping.ribo.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)
    
    @log
    def STAR(self):
        self.outPrefix = f'{self.outdir}/{self.sample}_'
        self.STAR_bam = f'{self.outdir}/{self.sample}_Aligned.sortedByCoord.out.bam'
        self.STAR_map_log = f'{self.outdir}/{self.sample}_Log.final.out'    
        cmd = ['STAR', '--runThreadN', str(self.thread), '--genomeDir', self.genomeDir,
            '--readFilesIn', self.fq, '--readFilesCommand', 'zcat', '--outFilterMultimapNmax',
            '1', '--outFileNamePrefix', self.outPrefix, '--outSAMtype', 'BAM', 'SortedByCoordinate',
            '--outFilterMatchNmin', self.outFilterMatchNmin]
        if self.out_unmapped:
            cmd += ['--outReadsUnmapped', 'Fastx']
        subprocess.check_call(cmd)

    @log
    def picard(self):
        self.picard_region_log = f'{self.outdir}/{self.sample}_region.log'
        cmd = [
            'picard',
            '-Xmx4G',
            '-XX:ParallelGCThreads=4',
            'CollectRnaSeqMetrics',
            'I=%s' %
            (self.STAR_bam),
            'O=%s' %
            (self.picard_region_log),
            'REF_FLAT=%s' %
            (self.refFlat),
            'STRAND=NONE',
            'VALIDATION_STRINGENCY=SILENT']
        res = subprocess.run(cmd, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)

    @log
    def run(self):
        self.STAR()
        self.picard()
        if self.debug:
            self.ribo()
        self.format_stat()
        self.report()
    
    def report(self):
        t = reporter(
            name='STAR',
            assay=self.assay,
            sample=self.sample,
            stat_file=self.stats_file,
            outdir=self.outdir + '/..',
            plot=self.plot)
        t.get_report()


def STAR(args):
    mapping = Step_mapping(
        args.sample, 
        args.outdir, 
        args.assay, 
        args.thread,
        args.fq, 
        args.genomeDir, 
        args.out_unmapped, 
        args.debug,
        args.outFilterMatchNmin)
    mapping.run()


def get_opts_STAR(parser, sub_program):
    if sub_program:
        parser.add_argument('--fq', required=True)
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument('--thread', default=1)
        parser.add_argument('--assay', help='assay', required=True)
        parser.add_argument('--debug', help='debug', action='store_true')
    parser.add_argument('--outFilterMatchNmin', help='STAR outFilterMatchNmin', default=0)
    parser.add_argument(
        '--out_unmapped',
        help='out_unmapped',
        action='store_true')
    parser.add_argument('--genomeDir', help='genome directory', required=True)
