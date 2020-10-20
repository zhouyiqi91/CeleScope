import argparse
import os
from celescope.tools.utils import parse_match_dir

CITESEQ_DIR = os.path.dirname(__file__)


def get_opts_analysis_cite(parser, sub_program):
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument('--match_dir', help='match_dir', required=True)
        parser.add_argument('--citeseq_mtx', help='citeseq matrix .gz file', required=True)
        parser.add_argument('--assay', help='assay', required=True)


def analysis_cite(args):

    if not os.path.exists(args.outdir):
        os.system('mkdir -p %s' % args.outdir)
        
    rds = parse_match_dir(args.match_dir)['rds']
    app = CITESEQ_DIR + "/analysis_cite.R"
    cmd = (
        f'Rscript {app} '
        f'--rds {rds} '
        f'--citeseq_mtx {args.citeseq_mtx} '
        f'--outdir {args.outdir} '
        f'--sample {args.sample} '
    )
    os.system(cmd)
    