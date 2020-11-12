#!/bin/env python
# coding=utf8

import os
import pandas as pd
import sys
import logging
from celescope.__init__ import __VERSION__, ASSAY_DICT
from celescope.tools.utils import log
from celescope.tools.report import reporter
from celescope.tools.__init__ import __PATTERN_DICT__
from .Chemistry import Chemistry


@log
def sample_info(args):

    sample = args.sample
    ASSAY = ASSAY_DICT[args.assay]
    version = __VERSION__
    outdir = args.outdir
    chemistry = args.chemistry
    fq1 = args.fq1
    fq1_file0 = fq1.split(',')[0]
    #transcriptome = args.genomeDir.split("/")[-1]

    # get chemistry
    if chemistry == 'auto':
        ch = Chemistry(fq1_file0)
        chemistry = ch.get_chemistry()

    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % outdir)

    stat = pd.DataFrame({
        "item": ["Sample ID", "Assay", "Chemistry", "Software Version"],
        "count": [sample, ASSAY, chemistry, version],
        },
        columns=["item", "count"]
    )
    stat_file = outdir + "/stat.txt"
    stat.to_csv(stat_file, sep=":", header=None, index=False)

    t = reporter(
        name='sample',
        assay=args.assay,
        sample=args.sample,
        stat_file=stat_file,
        outdir=outdir + '/..')
    t.get_report()
    return chemistry


def get_opts_sample(parser, sub_program):
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument('--assay', help='assay', required=True)
        parser.add_argument('--fq1', help='read1 fq file', required=True)
    parser.add_argument('--chemistry', choices=__PATTERN_DICT__.keys(), help='chemistry version', default='auto')
    
