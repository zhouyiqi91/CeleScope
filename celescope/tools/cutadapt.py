#!/bin/env python
# coding=utf8

import os
import sys
import re
import subprocess
import logging
from itertools import islice
import pandas as pd

from celescope.tools.utils import format_number, log
from celescope.tools.report import reporter


def format_stat(cutadapt_log, samplename):
    fh = open(cutadapt_log, 'r')
    stat_file = os.path.dirname(cutadapt_log) + '/stat.txt'
    # Total reads processed:...Total written (filtered):
    content = islice(fh, 9, 16)
    p_list = []
    for line in content:
        if line.strip() == '':
            continue
        line = re.sub(r'\s{2,}', r'', line)
        line = re.sub(r' bp', r'', line)
        line = re.sub(r'(?<=\d)\s+\(', r'(', line)
        line = line.strip()
        attr = line.split(":")
        p_list.append({"item": attr[0], "value": attr[1]})
    p_df = pd.DataFrame(p_list)
    p_df.iloc[0, 0] = 'Reads with Adapters'
    p_df.iloc[1, 0] = 'Reads too Short'
    p_df.iloc[2, 0] = 'Reads Written'
    p_df.iloc[3, 0] = 'Base Pairs Processed'
    p_df.iloc[4, 0] = 'Base Pairs Quality-Trimmed'
    p_df.iloc[5, 0] = 'Base Pairs Written'
    p_df.to_csv(stat_file, sep=':', index=False, header=None)

    fh.close()


@log
def cutadapt(args):
    # check dir
    if not os.path.exists(args.outdir):
        os.system('mkdir -p %s' % (args.outdir))

    # run cutadapt
    adapt = []
    for a in args.adapt:
        adapt.append('-a')
        adapt.append(a)

    out_fq2 = args.outdir + '/' + args.sample + '_clean_2.fq.gz'
    cmd = ['cutadapt'] + adapt + ['-n',
                                  str(len(args.adapt)),
                                  '-j',
                                  str(args.thread),
                                  '-m',
                                  str(args.minimum_length),
                                  '--nextseq-trim=' + str(args.nextseq_trim),
                                  '--overlap',
                                  str(args.overlap),
                                  '-o',
                                  out_fq2,
                                  args.fq]
    cutadapt.logger.info('%s' % (' '.join(cmd)))
    res = subprocess.run(cmd, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
    with open(args.outdir + '/cutadapt.log', 'wb') as fh:
        fh.write(res.stdout)

    format_stat(args.outdir + '/cutadapt.log', args.sample)
    t = reporter(
        name='cutadapt',
        assay=args.assay,
        sample=args.sample,
        stat_file=args.outdir +
        '/stat.txt',
        outdir=args.outdir +
        '/..')
    t.get_report()


def get_opts_cutadapt(parser, sub_program):
    if sub_program:
        parser.add_argument('--fq', help='fq file', required=True)
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument('--assay', help='assay', required=True)
    parser.add_argument(
        '--adapt',
        action='append',
        default=[
            'polyT=A{18}',
            'p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
        ])
    parser.add_argument(
        '--minimum_length',
        dest='minimum_length',
        help='minimum_length, default=20',
        default=20)
    parser.add_argument(
        '--nextseq-trim',
        dest='nextseq_trim',
        help='nextseq_trim, default=20',
        default=20)
    parser.add_argument(
        '--overlap',
        help='minimum overlap length, default=5',
        default=5)
    parser.add_argument('--thread', default=2)
