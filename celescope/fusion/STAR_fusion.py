import os
from celescope.tools.utils import log, STAR_util


@log
def STAR_fusion(args):
    STAR_util(
        args.sample,
        args.outdir,
        args.input_read,
        args.genomeDir,
        args.thread,
        args.outFilterMatchNmin,
    )


def get_opts_STAR_fusion(parser, sub_program):
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument("--input_read", required=True)
        parser.add_argument('--assay', help='assay', required=True)
    parser.add_argument(
        '--genomeDir',
        help='fusion gene STAR index genome directory',
        required=True)
    parser.add_argument("--thread", help='STAR thread', default=1)
    parser.add_argument("--outFilterMatchNmin", help='STAR outFilterMatchNmin', default=35)
