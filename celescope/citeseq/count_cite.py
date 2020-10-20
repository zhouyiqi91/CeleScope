import argparse
from .Count_cite import Count_cite

def get_opts_count_cite(parser, sub_program):
    parser.add_argument("--match_dir", help="matched scRNA-Seq CeleScope directory path", required=True)
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument('--assay', help='assay', required=True)
        parser.add_argument("--read_count_file", help="tag read count file")


def count_cite(args):

    count_cite_object = Count_cite(
        args.sample,
        args.outdir,
        args.assay,
        args.read_count_file,
        args.match_dir,
    )
    count_cite_object.run()
    count_cite_object.report()