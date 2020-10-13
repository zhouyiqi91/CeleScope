import logging
import os
import glob
import sys
import re
import gzip
import pandas as pd
import numpy as np
import subprocess
import time
import argparse
import pysam
from datetime import timedelta
from collections import defaultdict
from functools import wraps
from collections import Counter
import json
import celescope.tools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from celescope.tools.__init__ import __PATTERN_DICT__

tools_dir = os.path.dirname(celescope.tools.__file__)


def log(func):
    '''
    return logger.
    logging start and done.
    '''
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    module = func.__module__
    name = func.__name__
    logger_name = f'{module}.{name}'
    logger = logging.getLogger(logger_name)

    @wraps(func)
    def wrapper(*args, **kwargs):
        logger.info('start...')
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        used = timedelta(seconds=end - start)
        logger.info(f'done. time used: {used}')
        return result

    wrapper.logger = logger
    return wrapper


def arg_str(arg, arg_name):
    '''
    return action store_true arguments as string
    '''
    if arg:
        return '--' + arg_name
    return ''


def read_barcode_file(match_dir):
    '''
    multi version compatible
    '''
    match_barcode_file1 = glob.glob(
        f"{match_dir}/*count/*_cellbarcode.tsv")
    match_barcode_file2 = glob.glob(
        f"{match_dir}/*count/*matrix_10X/*_cellbarcode.tsv")
    match_barcode_file3 = glob.glob(
        f"{match_dir}/*count/*matrix_10X/*barcodes.tsv")
    match_barcode_file = (
        match_barcode_file1 +
        match_barcode_file2 +
        match_barcode_file3)[0]
    match_barcode, cell_total = read_one_col(match_barcode_file)
    match_barcode = set(match_barcode)
    return match_barcode, cell_total


def format_stat(count, total_count):
    percent = round(count / total_count * 100, 2)
    string = f'{count}({percent}%)'
    return string


def multi_opts(assay):
    readme = f'{assay} multi-samples'
    parser = argparse.ArgumentParser(readme)
    parser.add_argument('--mod', help='mod, sjm or shell', choices=['sjm', 'shell'], default='sjm')
    parser.add_argument(
        '--mapfile',
        help='''
            tsv file, 4 columns:
            1st col: LibName;
            2nd col: DataDir;
            3rd col: SampleName;
            4th col: Cell number or match_dir, optional;
        ''',
        required=True)
    parser.add_argument('--chemistry', choices=__PATTERN_DICT__.keys(), help='chemistry version')
    parser.add_argument('--whitelist', help='cellbarcode list')
    parser.add_argument('--linker', help='linker')
    parser.add_argument('--pattern', help='read1 pattern')
    parser.add_argument('--outdir', help='output dir', default="./")
    parser.add_argument(
        '--adapt',
        action='append',
        help='adapter sequence',
        default=[
            'polyT=A{15}',
            'p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'])
    parser.add_argument(
        '--minimum-length',
        dest='minimum_length',
        help='minimum_length',
        default=20)
    parser.add_argument(
        '--nextseq-trim',
        dest='nextseq_trim',
        help='nextseq_trim',
        default=20)
    parser.add_argument(
        '--overlap',
        help='minimum overlap length, default=5',
        default=5)
    parser.add_argument(
        '--lowQual',
        type=int,
        help='max phred of base as lowQual',
        default=0)
    parser.add_argument(
        '--lowNum',
        type=int,
        help='max number with lowQual allowed',
        default=2)
    parser.add_argument(
        '--rm_files',
        action='store_true',
        help='remove redundant fq.gz and bam after running')
    parser.add_argument('--steps_run', help='steps to run', default='all')
    return parser


def link_data(outdir, fq_dict):
    raw_dir = f'{outdir}/data_give/rawdata'
    os.system('mkdir -p %s' % (raw_dir))
    with open(raw_dir + '/ln.sh', 'w') as fh:
        fh.write('cd %s\n' % (raw_dir))
        for s, arr in fq_dict.items():
            fh.write('ln -sf %s %s\n' % (arr[0], s + '_1.fq.gz'))
            fh.write('ln -sf %s %s\n' % (arr[1], s + '_2.fq.gz'))


def gene_convert(gtf_file):

    gene_id_pattern = re.compile(r'gene_id "(\S+)";')
    gene_name_pattern = re.compile(r'gene_name "(\S+)"')
    id_name = {}
    c = Counter()
    with open(gtf_file) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('#'):
                continue
            tabs = line.split('\t')
            gtf_type, attributes = tabs[2], tabs[-1]
            if gtf_type == 'gene':
                gene_id = gene_id_pattern.findall(attributes)[-1]
                gene_name = gene_name_pattern.findall(attributes)[-1]
                c[gene_name] += 1
                if c[gene_name] > 1:
                    gene_name = f'{gene_name}_{c[gene_name]}'
                id_name[gene_id] = gene_name
    return id_name


@log
def process_read(
    read2_file, pattern_dict, barcode_dict, linker_dict,
    barcode_length, linker_length):

    # if valid, return (True)
    metrics = defaultdict(int)
    res_dict = genDict(dim=3)
    read2 = gzip.open(read2_file, "rt")
    while True:
        line1 = read2.readline()
        line2 = read2.readline()
        line3 = read2.readline()
        line4 = read2.readline()
        if not line4:
            break
        metrics['Total Reads'] += 1
        attr = str(line1).strip("@").split("_")
        barcode = str(attr[0])
        umi = str(attr[1])
        seq = line2.strip()
        if linker_dict:
            seq_linker = seq_ranges(seq, pattern_dict['L'])
            if len(seq_linker) < linker_length:
                metrics['Reads Unmapped too Short'] += 1
                continue
        if barcode_dict:
            seq_barcode = seq_ranges(seq, pattern_dict['C'])
            if barcode_length != len(seq_barcode):
                miss_length = barcode_length - len(seq_barcode)
                if miss_length > 2:
                    metrics['Reads Unmapped too Short'] += 1
                    continue
                seq_barcode = seq_barcode + "A" * miss_length                    
        
        # check linker
        valid_linker = False
        for linker_name in linker_dict:
            if hamming_correct(linker_dict[linker_name], seq_linker):
                valid_linker = True
                break

        if not valid_linker:
            metrics['Reads Unmapped Invalid Linker'] += 1
            continue

        # check barcode
        valid_barcode = False
        for barcode_name in barcode_dict:
            if hamming_correct(barcode_dict[barcode_name], seq_barcode):
                res_dict[barcode][barcode_name][umi] += 1
                valid_barcode = True
                break

        if not valid_barcode:
            metrics['Reads Unmapped Invalid Barcode'] += 1
            continue

        # mapped
        metrics['Reads Mapped'] += 1
        if metrics['Reads Mapped'] % 1000000 == 0:
            process_read.logger.info(str(metrics['Reads Mapped']) + " reads done.")

    return res_dict, metrics


def seq_ranges_exception(seq, pattern_dict):
    # get subseq with intervals in arr and concatenate
    length = len(seq)
    for x in pattern_dict:
        if length < x[1]:
            raise Exception(f'invalid seq range {x[0]}:{x[1]} in read')
    return ''.join([seq[x[0]:x[1]]for x in pattern_dict])


def seq_ranges(seq, pattern_dict):
    # get subseq with intervals in arr and concatenate
    return ''.join([seq[x[0]:x[1]]for x in pattern_dict])


def read_one_col(file):
    df = pd.read_csv(file, header=None)
    col1 = list(df.iloc[:, 0])
    num = len(col1)
    return col1, num


def read_fasta(fasta_file, equal=False):
    # seq must have equal length
    dict = {}
    LENGTH = 0
    with open(fasta_file, "rt") as f:
        while True:
            line1 = f.readline()
            line2 = f.readline()
            if not line2:
                break
            name = line1.strip(">").strip()
            seq = line2.strip()
            seq_length = len(seq)
            if equal:
                if LENGTH == 0:
                    LENGTH = seq_length
                else:
                    if LENGTH != seq_length:
                        raise Exception(f"{fasta_file} have different seq length")
            dict[name] = seq
    if equal:
        return dict, LENGTH
    return dict


def hamming_correct(string1, string2):
    threshold = len(string1) / 10 + 1
    if hamming_distance(string1, string2) < threshold:
        return True
    return False


def hamming_distance(string1, string2):
    distance = 0
    length = len(string1)
    length2 = len(string2)
    if (length != length2):
        raise Exception("string1 and string2 do not have same length")
    for i in range(length):
        if string1[i] != string2[i]:
            distance += 1
    return distance


def gen_stat(df, stat_file):
    # 3cols: item count total_count

    def add_percent(row):
        count = row['count']
        percent = count / row['total_count']
        value = f'{format_number(count)}({round(percent * 100, 2)}%)'
        return value
    df.loc[~df['total_count'].isna(), 'value'] = df.loc[~df['total_count'].isna(), :].apply(
        lambda row: add_percent(row), axis=1
    )
    df.loc[df['total_count'].isna(), 'value'] = df.loc[df['total_count'].isna(), :].apply(
        lambda row: f'{format_number(row["count"])}', axis=1
    )
    df = df.loc[:, ["item", "value"]]
    df.to_csv(stat_file, sep=":", header=None, index=False)


def get_fq(library_id, library_path):
    try:
        pattern1_1 = library_path + '/' + library_id + '*' + '_1.fq.gz'
        pattern1_2 = f'{library_path}/*{library_id}*R1.fastq.gz'
        pattern1_3 = f'{library_path}/*{library_id}*R1_001.fastq.gz'
        pattern2_1 = library_path + '/' + library_id + '*' + '_2.fq.gz'
        pattern2_2 = f'{library_path}/*{library_id}*R2.fastq.gz'
        pattern2_3 = f'{library_path}/*{library_id}*R2_001.fastq.gz'
        fq1 = ",".join(glob.glob(pattern1_1) + glob.glob(pattern1_2) + glob.glob(pattern1_3))
        fq2 = ",".join(glob.glob(pattern2_1) + glob.glob(pattern2_2) + glob.glob(pattern2_3))
        if len(fq1) == 0:
            sys.exit('Invalid fastq name pattern!')
    except IndexError as e:
        sys.exit("Mapfile Error:" + str(e))
    return fq1, fq2


def parse_map_col4(mapfile, default_val):
    fq_dict = defaultdict(list)
    col4_dict = defaultdict(list)
    with open(mapfile) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'):
                continue
            tmp = line.split()
            library_id = tmp[0]
            library_path = tmp[1]
            sample_name = tmp[2]
            if len(tmp) == 4:
                col4 = tmp[3]
            else:
                col4 = default_val
            fq1, fq2 = get_fq(library_id, library_path)

            if sample_name in fq_dict:
                fq_dict[sample_name][0].append(fq1)
                fq_dict[sample_name][1].append(fq2)
            else:
                fq_dict[sample_name] = [[fq1], [fq2]]
            col4_dict[sample_name] = col4

    for sample_name in fq_dict:
        fq_dict[sample_name][0] = ",".join(fq_dict[sample_name][0])
        fq_dict[sample_name][1] = ",".join(fq_dict[sample_name][1])

    return fq_dict, col4_dict


def generate_sjm(cmd, name, conda, m=1, x=1, steps_run='all'):
    res_cmd = f'''
job_begin
    name {name}
    sched_options -w n -cwd -V -l vf={m}g,p={x}
    cmd source activate {conda}; {cmd}
job_end
'''

    return res_cmd


def merge_report(
    fq_dict, steps, last_step, sjm_cmd,
    sjm_order, logdir, conda, outdir, rm_files):    
    step = "merge_report"
    steps_str = ",".join(steps)
    samples = ','.join(fq_dict.keys())
    app = tools_dir + '/merge_table.py'
    cmd = (
        f'python {app} --samples {samples} '
        f'--steps {steps_str} --outdir {outdir}'
    )
    if rm_files:
        cmd += ' --rm_files'
    sjm_cmd += generate_sjm(cmd, 'merge_report', conda)
    for sample in fq_dict:
        sjm_order += f'order {step} after {last_step}_{sample}\n'
    with open(logdir + '/sjm.job', 'w') as fh:
        fh.write(sjm_cmd + '\n')
        fh.write(sjm_order)


def format_number(number: int) -> str:
    return format(number, ",")


@log
def glob_genomeDir(genomeDir, fa=False):
    refFlat = glob.glob(genomeDir + "/*.refFlat")
    if (len(refFlat) > 1):
        sys.exit("ERROR: Multiple refFlat file in " + genomeDir)
    elif (len(refFlat) == 0):
        sys.exit("ERROR: refFlat file not found in " + genomeDir)
    else:
        refFlat = refFlat[0]
        glob_genomeDir.logger.info("refFlat file found: " + refFlat)

    gtf = glob.glob(genomeDir + "/*.gtf")
    if (len(gtf) == 0):
        sys.exit("ERROR: gtf file not found in " + genomeDir)
    elif (len(gtf) > 1):
        gtf = glob.glob(genomeDir + "/*.chr.gtf")
        if (len(gtf) == 0):
            sys.exit("ERROR: No chr gtf file in "+ genomeDir)
        if (len(gtf) > 1):
            sys.exit("ERROR: Multiple gtf file in " + genomeDir)
        else:
            gtf = gtf[0]
            glob_genomeDir.logger.info("chr gtf file found: " + gtf)
    else:
        gtf = gtf[0]
        glob_genomeDir.logger.info("gtf file found: " + gtf)
    
    if fa:
        fasta = glob.glob(genomeDir + "/*.fa") + glob.glob(genomeDir + "/*.fasta")
        if len(fasta) > 1:
            sys.exit("ERROR: Multiple fasta file in " + genomeDir)
        fasta = fasta[0]
        return refFlat, gtf, fasta
    return refFlat, gtf


def barcode_filter_with_magnitude(
        df, plot='magnitude.pdf', col='UMI', percent=0.1, expected_cell_num=3000):
    # col can be readcount or UMI
    # determine validated barcodes
    df = df.sort_values(col, ascending=False)
    idx = int(expected_cell_num * 0.01) - 1
    idx = max(0, idx)

    # calculate read counts threshold
    threshold = int(df.iloc[idx, df.columns == col] * 0.1)
    threshold = max(1, threshold)
    validated_barcodes = df[df[col] > threshold].index

    fig = plt.figure()
    plt.plot(df[col])
    plt.hlines(threshold, 0, len(validated_barcodes), linestyle='dashed')
    plt.vlines(len(validated_barcodes), 0, threshold, linestyle='dashed')
    plt.title('expected cell num: %s\n%s threshold: %s\ncell num: %s' %
              (expected_cell_num, col, threshold, len(validated_barcodes)))
    plt.loglog()
    plt.savefig(plot)

    return (validated_barcodes, threshold, len(validated_barcodes))


def barcode_filter_with_kde(df, plot='kde.pdf', col='UMI'):
    # col can be readcount or UMI
    # filter low values
    df = df.sort_values(col, ascending=False)
    arr = np.log10([i for i in df[col] if i / float(df[col][0]) > 0.001])

    # kde
    x_grid = np.linspace(min(arr), max(arr), 10000)
    density = gaussian_kde(arr, bw_method=0.1)
    y = density(x_grid)

    local_mins = argrelextrema(y, np.less)
    log_threshold = x_grid[local_mins[-1][0]]
    threshold = np.power(10, log_threshold)
    validated_barcodes = df[df[col] > threshold].index

    # plot
    fig, (ax1, ax2) = plt.subplots(2, figsize=(6.4, 10))
    ax1.plot(x_grid, y)
    #ax1.axhline(y[local_mins[-1][0]], -0.5, log_threshold, linestyle='dashed')
    ax1.vlines(log_threshold, 0, y[local_mins[-1][0]], linestyle='dashed')
    ax1.set_ylim(0, 0.3)

    ax2.plot(df[col])
    ax2.hlines(threshold, 0, len(validated_barcodes), linestyle='dashed')
    ax2.vlines(len(validated_barcodes), 0, threshold, linestyle='dashed')
    ax2.set_title('%s threshold: %s\ncell num: %s' %
                  (col, int(threshold), len(validated_barcodes)))
    ax2.loglog()
    plt.savefig(plot)

    return (validated_barcodes, threshold, len(validated_barcodes))


def get_slope(x, y, window=200, step=10):
    assert len(x) == len(y)
    start = 0
    last = len(x)
    res = [[], []]
    while True:
        end = start + window
        if end > last:
            break
        z = np.polyfit(x[start:end], y[start:end], 1)
        res[0].append(x[start])
        res[1].append(z[0])
        start += step
    return res


def barcode_filter_with_derivative(
        df, plot='derivative.pdf', col='UMI', window=500, step=5):
    # col can be readcount or UMI
    # filter low values
    df = df.sort_values(col, ascending=False)
    y = np.log10([i for i in df[col] if i / float(df[col][0]) > 0.001])
    x = np.log10(np.arange(len(y)) + 1)

    # derivative
    res = get_slope(x, y, window=window, step=step)
    res2 = get_slope(res[0], res[1], window=window, step=step)
    g0 = [res2[0][i] for i, j in enumerate(res2[1]) if j > 0]
    cell_num = int(np.power(10, g0[0]))
    threshold = df[col][cell_num]
    validated_barcodes = df.index[0:cell_num]

    # plot
    fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(6.4, 15))
    ax1.plot(res[0], res[1])

    ax2.plot(res2[0], res2[1])
    ax2.set_ylim(-1, 1)

    ax3.plot(df[col])
    ax3.hlines(threshold, 0, len(validated_barcodes), linestyle='dashed')
    ax3.vlines(len(validated_barcodes), 0, threshold, linestyle='dashed')
    ax3.set_title('%s threshold: %s\ncell num: %s' %
                  (col, int(threshold), len(validated_barcodes)))
    ax3.loglog()
    plt.savefig(plot)

    return (validated_barcodes, threshold, len(validated_barcodes))


def genDict(dim=3, valType=int):
    if dim == 1:
        return defaultdict(valType)
    else:
        return defaultdict(lambda: genDict(dim - 1, valType=valType))


@log
def downsample(bam, barcodes, percent):
    """
    calculate median_geneNum and saturation based on a given percent of reads

    Args:
        bam - bam file
        barcodes(set) - validated barcodes
        percent(float) - percent of reads in bam file.

    Returns:
        percent(float) - input percent
        median_geneNum(int) - median gene number
        saturation(float) - sequencing saturation

    Description:
        Sequencing Saturation = 1 - n_deduped_reads / n_reads.
        n_deduped_reads = Number of unique (valid cell-barcode, valid UMI, gene) combinations among confidently mapped reads.
        n_reads = Total number of confidently mapped, valid cell-barcode, valid UMI reads.
    """
    logging.info('working' + str(percent))
    cmd = ['samtools', 'view', '-s', str(percent), bam]
    p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    # nesting defaultdicts in an arbitrary depth
    readDict = genDict()

    n_reads = 0
    while True:
        line = p1.stdout.readline()
        if not line.strip():
            break
        tmp = line.strip().split()
        if tmp[-1].startswith('XT:Z:'):
            geneID = tmp[-1].replace('XT:Z:', '')
            cell_barcode, umi = tmp[0].split('_')[0:2]
            # filter invalid barcode
            if cell_barcode in barcodes:
                n_reads += 1
                readDict[cell_barcode][umi][geneID] += 1
    p1.stdout.close()

    geneNum_list = []
    n_deduped_reads = 0
    for cell_barcode in readDict:
        genes = set()
        for umi in readDict[cell_barcode]:
            for geneID in readDict[cell_barcode][umi]:
                genes.add(geneID)
                if readDict[cell_barcode][umi][geneID] == 1:
                    n_deduped_reads += 1
        geneNum_list.append(len(genes))

    median_geneNum = np.median(geneNum_list) if geneNum_list else 0
    saturation = (1 - float(n_deduped_reads) / n_reads) * 100

    return "%.2f\t%.2f\t%.2f\n" % (
        percent, median_geneNum, saturation), saturation


def cluster_tsne_list(tsne_df):
    """
    tSNE_1	tSNE_2	cluster Gene_Counts
    return data list
    """
    sum_df = tsne_df.groupby(["cluster"]).agg("count").iloc[:, 0]
    percent_df = sum_df.transform(lambda x: round(x / sum(x) * 100, 2))
    res = []
    for cluster in sorted(tsne_df.cluster.unique()):
        sub_df = tsne_df[tsne_df.cluster == cluster]
        name = "cluster {cluster}({percent}%)".format(
            cluster=cluster, percent=percent_df[cluster])
        tSNE_1 = list(sub_df.tSNE_1)
        tSNE_2 = list(sub_df.tSNE_2)
        res.append({"name": name, "tSNE_1": tSNE_1, "tSNE_2": tSNE_2})
    return res


def marker_table(marker_df):
    """
    return html code
    """
    marker_df = marker_df.loc[:, ["cluster", "gene",
                                  "avg_logFC", "pct.1", "pct.2", "p_val_adj"]]
    marker_gene_table = marker_df.to_html(
        escape=False,
        index=False,
        table_id="marker_gene_table",
        justify="center")
    return marker_gene_table


def report_prepare(outdir, **kwargs):
    json_file = outdir + '/../.data.json'
    if not os.path.exists(json_file):
        data = {}
    else:
        fh = open(json_file)
        data = json.load(fh)
        fh.close()

    for key in kwargs:
        data[key] = kwargs[key]

    with open(json_file, 'w') as fh:
        json.dump(data, fh)


def parse_vcf(vcf_file, cols=['chrom', 'pos', 'alleles'], infos=['DP', 'GENE', 'CELL']):
    vcf = pysam.VariantFile(vcf_file)
    df = pd.DataFrame(columns=[col.upper() for col in cols] + infos + ['GT'])
    rec_dict = {}
    for rec in vcf.fetch():

        for col in cols:
            rec_dict[col.upper()] = getattr(rec, col)
            if col == 'alleles':
                rec_dict['ALLELES'] = '-'.join(rec_dict['ALLELES'])
                
        for info in infos:
            rec_dict[info] = rec.info[info]

        rec_dict['GT'] = [s['GT'] for s in rec.samples.values()][0]

        df = df.append(pd.Series(rec_dict),ignore_index=True)
    return df

