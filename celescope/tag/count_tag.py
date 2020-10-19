import argparse
from .Count_tag import Count_tag

def get_opts_count_tag(parser, sub_program):
    parser.add_argument(
        "--UMI_min",
        help="cells have tag_UMI>=UMI_min are considered as valid cell",
        default="auto")
    parser.add_argument("--dim", help="tag dimension", default=1)
    parser.add_argument(
        "--SNR_min",
        help="minimum signal to noise ratio",
        default="auto")
    parser.add_argument(
        "--combine_cluster",
        help="conbine cluster tsv file",
        default=None)
    parser.add_argument("--match_dir", help="matched scRNA-Seq CeleScope directory path", required=True)
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument('--assay', help='assay', required=True)
        parser.add_argument("--read_count_file", help="tag read count file")


def count_tag(args):

    count_tag_object = Count_tag(
        args.sample,
        args.outdir,
        args.assay,
        args.read_count_file,
        args.match_dir,
        args.UMI_min,
        args.SNR_min,
        args.combine_cluster,
        args.dim,
    )
    count_tag_object.run()
    count_tag_object.report()