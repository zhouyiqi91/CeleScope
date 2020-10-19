from .Mapping_tag import Mapping_tag
import argparse

def get_opts_mapping_tag(parser, sub_program):
    parser.add_argument("--fq_pattern", help="read2 fastq pattern")
    parser.add_argument("--linker_fasta", help="linker fasta")
    parser.add_argument("--barcode_fasta", help="barcode fasta")
    if sub_program:
        parser.add_argument("--fq", help="clean read2", required=True)
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument('--assay', help='assay', required=True)


def mapping_tag(args):
    sample = args.sample
    outdir = args.outdir
    assay = args.assay
    fq = args.fq
    fq_pattern = args.fq_pattern
    linker_fasta = args.linker_fasta
    barcode_fasta = args.barcode_fasta

    mapping_tag_object = Mapping_tag(
        sample,
        outdir,
        assay,
        fq,
        fq_pattern,
        linker_fasta,
        barcode_fasta,
    )
    mapping_tag_object.run()
