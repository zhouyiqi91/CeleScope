import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import celescope
import pysam
import numpy as np
import pandas as pd
import logging
from celescope.tools.utils import format_number, log, read_barcode_file
from celescope.tools.utils import format_stat
from celescope.tools.utils import read_one_col, gene_convert, glob_genomeDir
from celescope.tools.report import reporter


@log
def split_bam(bam, barcodes, outdir, sample, gene_id_name_dic, min_query_length):
    '''
    input:
        bam: bam from feauturCounts
        barcodes: cell barcodes, set
        gene_id_name_dic: id name dic
        min_query_length: minimum query length

    ouput:
        bam_dict: assign reads to cell barcodes and UMI
        count_dict: UMI counts per cell
        index: assign index(1-based) to cells
    '''

    # init
    count_dict = defaultdict(dict)
    bam_dict = defaultdict(dict)
    index_dict = defaultdict(dict)
    cells_dir = f'{outdir}/cells/'

    # read bam and split
    split_bam.logger.info('reading bam...')
    samfile = pysam.AlignmentFile(bam, "rb")
    header = samfile.header
    barcodes = list(barcodes)
    barcodes.sort()
    for read in samfile:
        attr = read.query_name.split('_')
        barcode = attr[0]
        umi = attr[1]
        if not read.has_tag('XT'):
            continue
        gene = read.get_tag('XT')
        query_length = read.infer_query_length()
        if (barcode in barcodes) and (gene in gene_id_name_dic) and (query_length >= min_query_length):
            gene_name = gene_id_name_dic[gene]
            read.set_tag(tag='GN', value=gene_name, value_type='Z')
            index = barcodes.index(barcode) + 1
            read.set_tag(tag='CL', value=f'CELL{index}', value_type='Z')
            # keep one read for each UMI
            if umi not in bam_dict[barcode]:
                bam_dict[barcode][umi] = read
            # count
            if gene_name not in count_dict[barcode]:
                count_dict[barcode][gene_name] = {}
            if umi in count_dict[barcode][gene_name]:
                count_dict[barcode][gene_name][umi] += 1
            else:
                count_dict[barcode][gene_name][umi] = 1

    split_bam.logger.info('writing cell bam...')
    # write new bam
    index = 0
    for barcode in barcodes:
        # init
        index += 1
        index_dict[index]['barcode'] = barcode
        index_dict[index]['valid'] = False

        # out bam
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
    df_count = pd.DataFrame(columns=['barcode', 'gene', 'UMI', 'read_count'])
    for barcode in count_dict:
        for gene in count_dict[barcode]:
            for umi in count_dict[barcode][gene]:
                read_count = count_dict[barcode][gene][umi]
                df_count = df_count.append({
                    'barcode': barcode,
                    'gene': gene,
                    'UMI': umi,
                    'read_count': read_count,
                }, ignore_index=True)
    count_file = f'{outdir}/{sample}_count.tsv'
    df_count.to_csv(count_file, sep='\t', index=False)

    return index_file, count_file


def call_snp(index, outdir, fasta):

    logging.info(f'Processing Cell {index}..')
    bam = f'{outdir}/cells/cell{index}/cell{index}.bam'
    # sort
    sorted_bam = f'{outdir}/cells/cell{index}/cell{index}_sorted.bam'
    cmd_sort = (
        f'samtools sort {bam} -o {sorted_bam}'
    )
    os.system(cmd_sort)

    # mpileup
    out_bcf = f'{outdir}/cells/cell{index}/cell{index}.bcf'
    cmd_mpileup = (
        f'samtools mpileup -g -f '
        f'{fasta} '
        f'{sorted_bam} 2>/dev/null | '
        f'bcftools call -mv -Ob '
        f'-o {out_bcf} '
        f'>/dev/null 2>&1 '
    )
    os.system(cmd_mpileup)

    # view
    out_vcf = f'{outdir}/cells/cell{index}/cell{index}.vcf'
    cmd_view = (
        f'bcftools view {out_bcf} > {out_vcf}'
    )
    os.system(cmd_view)


def read_index(index_file):
    df_index = pd.read_csv(index_file, sep='\t', index_col=0, dtype=object)
    df_valid = df_index[df_index['valid'] == 'True']
    return df_index, df_valid


@log
def call_all_snp(index_file, outdir, thread, fasta):
    all_res = []
    _df_index, df_valid = read_index(index_file)
    index_arg = df_valid.index
    outdir_arg = [outdir] * len(index_arg)
    fasta_arg = [fasta] * len(index_arg)
    with ProcessPoolExecutor(thread) as pool:
        for res in pool.map(call_snp, index_arg, outdir_arg, fasta_arg):
            all_res.append(res)


def process_vcf_header(line, sample):
    if line.startswith('##'):
        if line.find('samtools') != -1 or line.find('bcftools') != -1:
            return
        return line
    else:
        items = line.split('\t')
        items[-1] = sample
        new_line = '\t'.join(items) + '\n'
        return new_line


@log
def summary(index_file, count_file, outdir, sample):
    # init
    number = 0
    Number_of_Match_Cells_with_SNP = 0
    SNP_count_dict = defaultdict(int)
    coord_gene_dict = defaultdict(dict)

    # read index
    df_index, df_valid = read_index(index_file)

    # out vcf
    out_vcf = open(f'{outdir}/{sample}.vcf', 'wt')
    for index in df_valid.index:
        vcf_coords_dict = {}
        number += 1
        cell_vcf_file = f'{outdir}/cells/cell{index}/cell{index}.vcf'
        # vcf coords
        with open(cell_vcf_file, 'rt') as f:
            for line in f:
                if line.startswith("#"):
                    # add vcf and bam header
                    if number == 1:
                        new_line = process_vcf_header(line, sample)
                        if new_line:
                            out_vcf.write(new_line)
                    continue
                if line:
                    items = line.split('\t')
                    items[7] += f';CELL={index}'
                    new_line = '\t'.join(items)
                    out_vcf.write(new_line)
                    chrom = str(items[0])
                    pos = int(items[1])
                    if chrom not in vcf_coords_dict:
                        vcf_coords_dict[chrom] = set([pos])
                    else:
                        vcf_coords_dict[chrom].add(pos)
                    SNP_count_dict[index] += 1

        # add bam header
        if number == 1:
            cell_bam_file = f'{outdir}/cells/cell{index}/cell{index}_sorted.bam'
            cell_bam = pysam.AlignmentFile(cell_bam_file, "rb")
            header = cell_bam.header
            out_bam = pysam.AlignmentFile(f'{outdir}/{sample}.bam', "wb", header=header)

        # add bam
        if len(vcf_coords_dict) > 0:
            Number_of_Match_Cells_with_SNP += 1
            cell_bam_file = f'{outdir}/cells/cell{index}/cell{index}_sorted.bam'
            cell_bam = pysam.AlignmentFile(cell_bam_file, "rb")
            for read in cell_bam:
                bam_ref = str(read.reference_name)
                gene_name = read.get_tag('GN')
                aligned_pairs = read.get_aligned_pairs()
                align_dict = {}
                for pair in aligned_pairs:
                    ref_pos = pair[1]
                    read_pos = pair[0]
                    if ref_pos:
                        align_dict[ref_pos] = read_pos
                if bam_ref in vcf_coords_dict.keys():
                    read_flag = False
                    for pos in vcf_coords_dict[bam_ref]:
                        if pos in align_dict:
                            read_flag = True
                            coord_gene_dict[bam_ref][pos] = gene_name
                    if read_flag:
                        out_bam.write(read)

    out_vcf.close()
    out_bam.close()
    pysam.sort("-o", f'{outdir}/{sample}_sorted.bam', f'{outdir}/{sample}.bam')
    cmd = f'samtools index {outdir}/{sample}_sorted.bam'
    os.system(cmd)

    # annotate vcf
    anno_vcf = open(f'{outdir}/{sample}_anno.vcf', 'wt')
    with open(f'{outdir}/{sample}.vcf', 'rt') as vcf:
        for line in vcf:
            if line.startswith('#'):
                anno_vcf.write(line)
                continue
            items = line.split('\t')
            chrom = str(items[0])
            pos = int(items[1])
            gene_name = coord_gene_dict[chrom][pos]
            items[7] += f';GENE={gene_name}'
            new_line = '\t'.join(items)
            anno_vcf.write(new_line)
    anno_vcf.close()

    # rm
    #os.remove(f'{outdir}/{sample}.vcf')
    #os.remove(f'{outdir}/{sample}.bam')

    # stat
    stats = pd.Series()
    n_match_cell = len(df_index.index)

    df_count = pd.read_csv(count_file, sep='\t')
    df_count_read = df_count.groupby('barcode').agg({'read_count':sum})
    read_total = sum(df_count_read['read_count'])
    Mean_Reads_per_Cell = round((read_total / n_match_cell), 2)
    stats = stats.append(pd.Series(
        Mean_Reads_per_Cell,
        index=['Mean Reads per Cell']
    ))
    df_count_UMI = df_count.groupby('barcode').agg({'UMI':'count'})
    UMI_total = sum(df_count_UMI['UMI'])
    Mean_UMIs_per_Cell = round((UMI_total / n_match_cell), 2)
    stats = stats.append(pd.Series(
        Mean_UMIs_per_Cell,
        index=['Mean UMIs per Cell']
    ))

    stats = stats.append(pd.Series(
        format_stat(Number_of_Match_Cells_with_SNP, n_match_cell),
        index=['Number of Cells with Variants']
    ))

    SNP_counts = list(SNP_count_dict.values())
    Mean_SNP_per_Cell = round(np.mean(SNP_counts), 3)
    stats = stats.append(pd.Series(
        Mean_SNP_per_Cell,
        index=['Mean Variants per Cell with Variants']
    ))

    stat_file = f'{outdir}/stat.txt'
    stats.to_csv(stat_file, sep=':', header=False)

    t = reporter(
        name='snpCalling', assay='snp', sample=sample,
        stat_file=stat_file, outdir=outdir + '/..')
    t.get_report()


@log
def convert(gene_list_file, gtf):
    gene_list_name, _count = read_one_col(gene_list_file)
    id_name = gene_convert(gtf)
    name_id = {}
    for id in id_name:
        name = id_name[id]
        name_id[name] = id
    gene_id_name_dic = {}
    for gene_name in gene_list_name:
        gene_id = name_id[gene_name]
        gene_id_name_dic[gene_id] = gene_name
    return gene_id_name_dic


@log
def snpCalling(args):

    sample = args.sample
    outdir = args.outdir
    thread = int(args.thread)
    match_dir = args.match_dir
    bam = args.bam
    genomeDir = args.genomeDir
    gene_list_file = args.gene_list
    min_query_length = args.min_query_length

    # process args
    barcodes, _nCell = read_barcode_file(match_dir)

    # check dir
    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % (outdir))

    # get genome file
    _refFlat, gtf, fasta = glob_genomeDir(genomeDir, fa=True)

    # convert gene
    gene_id_name_dic = convert(gene_list_file, gtf)

    # split bam
    index_file, count_file = split_bam(
        bam, barcodes, outdir, sample, gene_id_name_dic, min_query_length)

    # snp
    call_all_snp(index_file, outdir, thread, fasta)

    # summary
    summary(index_file, count_file, outdir, sample)


def get_opts_snpCalling(parser, sub_program):
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument('--assay', help='assay', required=True)
        parser.add_argument("--thread", help='number of thread', default=1)
        parser.add_argument("--bam", help='featureCounts bam', required=True)
        parser.add_argument("--genomeDir", help='genomeDir', required=True)
    parser.add_argument(
        "--match_dir", help="match scRNA-Seq dir", required=True)
    parser.add_argument("--gene_list", help='gene_list', required=True)
    parser.add_argument("--min_query_length", help='minimum query length', default=35)
