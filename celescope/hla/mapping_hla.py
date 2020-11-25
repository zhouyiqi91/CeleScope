import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import celescope
import pysam
import pandas as pd
from celescope.tools.utils import format_number, log, read_barcode_file


@log
def split_bam(bam,barcodes, outdir, sample):
    
    # init
    count_dict = defaultdict(dict)
    bam_dict = defaultdict(dict)
    index_dict = defaultdict(dict)
    cells_dir = f'{outdir}/cells/'

    # read bam and split
    split_bam.logger.info('reading bam...')
    samfile = pysam.AlignmentFile(bam, "rb")
    header = samfile.header
    for read in samfile:
        attr = read.query_name.split('_')
        barcode = attr[0]
        umi = attr[1]
        if barcode in barcodes:
            # keep one read for each UMI
            if umi not in bam_dict[barcode]:
                bam_dict[barcode][umi] = read
            if umi in count_dict[barcode]:
                count_dict[barcode][umi] += 1
            else:
                count_dict[barcode][umi] = 1

    split_bam.logger.info('writing cell bam...')
    # write new bam
    index = 0
    for barcode in barcodes:
        # init
        index += 1
        index_dict[index]['barcode'] = barcode
        index_dict[index]['valid'] = False

        # out
        if barcode in bam_dict:
            cell_dir = f'{cells_dir}/cell{index}'
            cell_bam_file = f'{cell_dir}/cell{index}.bam'
            if not os.path.exists(cell_dir):
                os.makedirs(cell_dir)
            index_dict[index]['valid'] = True
            cell_bam = pysam.AlignmentFile(
                f'{cell_bam_file}', "wb", header=header)
            for umi in bam_dict[barcode]:
                read = bam_dict[barcode][umi]
                cell_bam.write(read)
            cell_bam.close()

    # out df_index
    df_index = pd.DataFrame(index_dict).T
    df_index.index.name = 'cell_index'
    index_file = f'{outdir}/{sample}_cell_index.tsv'
    df_index.to_csv(index_file, sep='\t')

    # out count_dict
    df_count = pd.DataFrame(count_dict).T
    df_temp = df_count.stack()
    df_temp = df_temp.astype(int)
    df_temp = df_temp.reset_index()
    df_temp.columns = ['barcode', 'UMI', 'read_count']
    count_file = f'{outdir}/{sample}_UMI_count.tsv'
    df_temp.to_csv(count_file, sep='\t', index=False)

    return index_file, count_file
    


def sub_typing(bam):
    
    outdir = os.path.dirname(bam)
    cmd = (
        f'arcasHLA extract '
        f'{bam} '
        f'--outdir {outdir} '
    )
    os.system(cmd)

def sub_typing2(fastq):
    
    outdir = os.path.join(os.path.dirname(fastq), "../..")
    
    cmd = (
        f'arcasHLA genotype '
        f'--min_count 1 '
        f'--drop_iterations 1 '
        f'{fastq} '  
        f'--outdir {outdir}/genotype/ '
    )
    os.system(cmd)
    


def read_index(index_file):
    df_index = pd.read_csv(index_file, sep='\t', index_col=0, dtype=object)
    df_valid = df_index[df_index['valid'] == 'True']
    return df_valid


@log
def hla_typing(index_file, outdir, thread):
    all_res = []
    df_valid = read_index(index_file)

    bams = [
        f'{outdir}/cells/cell{index}/cell{index}.bam' for index in df_valid.index]
    
    with ProcessPoolExecutor(thread) as pool:
        for res in pool.map(sub_typing,bams):
            all_res.append(res)

@log
def hla_typing2(index_file, outdir, thread):
    all_res = []
    df_valid = read_index(index_file)
    out_type_dir = f'{outdir}/genotype/'
    if not os.path.exists(out_type_dir):
        os.mkdir(out_type_dir)
    
    fastqs = [
        f'{outdir}/cells/cell{index}/cell{index}.extracted.fq.gz' for index in df_valid.index]
    with ProcessPoolExecutor(thread) as pool:
        for res in pool.map(sub_typing2,fastqs):
            all_res.append(res)

@log
def summary(outdir):
    
    cmd = (
        f'arcasHLA merge '
        f'-i {outdir}/genotype/ '
        f'-o {outdir} '
    )
    os.system(cmd)


@log
def mapping_hla(args):

    sample = args.sample
    outdir = args.outdir
    bam = args.bam
    thread = int(args.thread)
    match_dir = args.match_dir

    # process args
    barcodes, _nCell = read_barcode_file(match_dir)

    # check dir
    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % (outdir))

   

    # split bam
    index_file, _count_file = split_bam(bam,barcodes, outdir, sample)

    # typing
    hla_typing(index_file, outdir, thread)

    hla_typing2(index_file, outdir, thread)

    # summary
    summary(outdir)


def get_opts_mapping_hla(parser, sub_program):
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument("--bam", required=True)
        parser.add_argument('--assay', help='assay', required=True)
    parser.add_argument(
        "--match_dir", help="match scRNA-Seq dir", required=True)
    parser.add_argument("--thread", help='number of thread', default=1)
