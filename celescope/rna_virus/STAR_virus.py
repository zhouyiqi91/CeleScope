#!/bin/env python
# coding=utf8

import os
import logging
logger1 = logging.getLogger(__name__)


def STAR_virus(args):

    logger1.info('align virus genome...')

    sample = args.sample
    outdir = args.outdir
    input_read = args.input_read
    virus_genomeDir = args.virus_genomeDir
    thread = args.thread
    outFilterMatchNmin = args.outFilterMatchNmin

    # check dir
    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % (outdir))

    out_prefix = outdir + "/" + sample + "_virus_"
    out_BAM = out_prefix + "Aligned.sortedByCoord.out.bam"

    # host genome align
    cmd = "STAR \
 --genomeDir {genome} \
 --readFilesIn {input_read}\
 --outFilterMatchNmin 35\
 --outSAMtype BAM SortedByCoordinate\
 --runThreadN {runThreadN}\
 --limitBAMsortRAM 1933647594\
 --outFileNamePrefix {out_prefix}".format(genome=virus_genomeDir,
                                          input_read=input_read, runThreadN=thread, out_prefix=out_prefix)

    # add gz
    if input_read[len(input_read) - 2:] == "gz":
        cmd += ' --readFilesCommand zcat'

    logger1.info(cmd)
    os.system(cmd)
    logger1.info("align virus genome done.")

    cmd = "samtools index {out_BAM}".format(out_BAM=out_BAM)
    logger1.info(cmd)
    os.system(cmd)
    logger1.info("index done.")


def get_opts_STAR_virus(parser, sub_program):
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument("--input_read", required=True)
        parser.add_argument("--thread", help='STAR thread', default=1)
        parser.add_argument('--assay', help='assay', required=True)
    parser.add_argument(
        '--virus_genomeDir',
        help='virus genome dir',
        required=True)
    parser.add_argument("--outFilterMatchNmin", help='STAR outFilterMatchNmin', default=35)
