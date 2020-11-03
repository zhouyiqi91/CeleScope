#!/bin/env python
# coding=utf8
import os
import re
import io
import gzip
import subprocess
import sys
import glob
import pandas as pd
from collections import defaultdict, Counter
from itertools import combinations, permutations, islice
from xopen import xopen
from celescope.tools.utils import format_number, log, seq_ranges, read_fasta, genDict
from celescope.tools.report import reporter
from celescope.tools.__init__ import __PATTERN_DICT__

barcode_corrected_num = 0

# 定义输出格式
stat_info = '''
    Raw Reads: %s
    Valid Reads: %s(%s)
    Q30 of Barcodes: %.2f%%
    Q30 of UMIs: %.2f%%
'''


def ord2chr(q, offset=33):
    return chr(int(q) + offset)


# 生成错配字典
def generate_mis_seq(seq, n=1, bases='ACGTN'):
    # 以随机bases中的碱基替换seq中的n个位置，产生的错配字典
    # 返回字典，错配序列为key，
    # (正确序列，错配碱基数目，错配碱基位置，原始碱基，新碱基)组成的元组
    # 作为字典的值

    length = len(seq)
    assert length >= n, "err number should not be larger than sequence length!"
    res = {}
    seq_arr = list(seq)
    pos_group = list(combinations(range(0, length), n))
    bases_group = list(permutations(bases, n))

    for g in pos_group:
        for b in bases_group:
            seq_tmp = seq_arr[:]
            mis_num = n
            raw_tmp = []
            for i in range(n):
                raw_base = seq_tmp[g[i]]
                new_base = b[i]

                if raw_base == new_base:
                    mis_num -= 1

                raw_tmp.append(raw_base)
                seq_tmp[g[i]] = new_base

            if mis_num != 0:
                res[''.join(seq_tmp)] = (seq, mis_num, ','.join(
                    [str(i) for i in g]), ','.join(raw_tmp), ','.join(b))
    return(res)


@log
def generate_seq_dict(seqlist, n=1):
    seq_dict = {}
    with open(seqlist, 'r') as fh:
        for seq in fh:
            seq = seq.strip()
            if seq == '':
                continue
            seq_dict[seq] = (seq, 0, -1, 'X', 'X')
            for k, v in generate_mis_seq(seq, n).items():
                # duplicate key
                if k in seq_dict:
                    generate_seq_dict.logger.warning('barcode %s, %s\n%s, %s' %
                                                     (v, k, seq_dict[k], k))
                else:
                    seq_dict[k] = v
    return seq_dict


@log
def parse_pattern(pattern):
    # 解析接头结构，返回接头结构字典
    # key: 字母表示的接头, value: 碱基区间列表
    # eg.: C8L10C8L10C8U8T30
    # defaultdict(<type 'list'>:
    # {'C': [[0, 8], [18, 26], [36, 44]], 'U': [[44, 52]], 'L': [[8, 18], [26, 36]], 'T': [[52, 82]]})
    pattern_dict = defaultdict(list)
    p = re.compile(r'([CLUNT])(\d+)')
    tmp = p.findall(pattern)
    if not tmp:
        parse_pattern.logger.error(f'Invalid pattern: {pattern}')
        sys.exit()
    start = 0
    for item in tmp:
        end = start + int(item[1])
        pattern_dict[item[0]].append([start, end])
        start = end
    return pattern_dict


def get_scope_bc(bctype):
    import celescope
    root_path = os.path.dirname(celescope.__file__)
    linker_f = glob.glob(f'{root_path}/data/{bctype}/linker*')[0]
    whitelist_f = f'{root_path}/data/{bctype}/bclist'
    return linker_f, whitelist_f


def read_fastq(f):
    """
    Return tuples: (name, sequence, qualities).
    qualities is a string and it contains the unmodified, encoded qualities.
    """
    i = 3
    for i, line in enumerate(f):
        if i % 4 == 0:
            assert line.startswith('@'), ("Line {0} in FASTQ file is expected to start with '@', "
                                          "but found {1!r}".format(i + 1, line[:10]))
            name = line.strip()[1:]
        elif i % 4 == 1:
            sequence = line.strip()
        elif i % 4 == 2:
            line = line.strip()
            assert line.startswith('+'), ("Line {0} in FASTQ file is expected to start with '+', "
                                          "but found {1!r}".format(i + 1, line[:10]))
        elif i % 4 == 3:
            qualities = line.rstrip('\n\r')
            yield name, sequence, qualities
    if i % 4 != 3:
        raise Exception("FASTQ file ended prematurely")


def low_qual(quals, minQ='/', num=2):
    # print(ord('/')-33)           14
    return True if len([q for q in quals if q < minQ]) > num else False


def no_polyT(seq, strictT=0, minT=10):
    # strictT requires the first nth position to be T
    if seq[:strictT] != 'T' * strictT or seq.count('T') < minT:
        return True
    else:
        return False


def no_barcode(seq_arr, mis_dict, err_tolerance=1):
    global barcode_corrected_num
    tmp = [mis_dict[seq][0:2] if seq in mis_dict else (
        'X', 100) for seq in seq_arr]
    err = sum([t[1] for t in tmp])
    if err > err_tolerance:
        return True
    else:
        if err > 0:
            barcode_corrected_num += 1
            return ''.join([t[0] for t in tmp])
        else:
            return "correct"


def no_linker(seq, linker_dict):
    return False if seq in linker_dict else True


def check_seq(seq_file, pattern_dict, seq_abbr):
    length = 0
    for item in pattern_dict[seq_abbr]:
        start = item[0]
        end = item[1]
        length += end - start
    with open(seq_file, 'r') as fh:
        for seq in fh:
            seq = seq.strip()
            if seq == '':
                continue
            if len(seq) != length:
                raise Exception(
                    f'length of L in pattern ({length}) do not equal to length in {seq_file} ({len(seq)}) !')


@log
def barcode(args):

    # check dir
    if not os.path.exists(args.outdir):
        os.system('mkdir -p %s' % args.outdir)

    bc_pattern = __PATTERN_DICT__[args.chemistry]
    if (bc_pattern):
        (linker, whitelist) = get_scope_bc(args.chemistry)
    else:
        bc_pattern = args.pattern
        linker = args.linker
        whitelist = args.whitelist
    if (not linker) or (not whitelist) or (not bc_pattern):
        barcode.logger.error("invalid chemistry or [pattern,linker,whitelist]")
        sys.exit()

    # parse pattern to dict, C8L10C8L10C8U8
    # defaultdict(<type 'list'>, {'C': [[0, 8], [18, 26], [36, 44]], 'U':
    # [[44, 52]], 'L': [[8, 18], [26, 36]]})
    pattern_dict = parse_pattern(bc_pattern)

    # check linker
    check_seq(linker, pattern_dict, "L")

    bool_T = True if 'T' in pattern_dict else False
    bool_L = True if 'L' in pattern_dict else False

    C_len = sum([item[1] - item[0] for item in pattern_dict['C']])

    barcode_qual_Counter = Counter()
    umi_qual_Counter = Counter()
    C_U_base_Counter = Counter()
    args.lowQual = ord2chr(args.lowQual)

    # generate list with mismatch 1, substitute one base in raw sequence with
    # A,T,C,G
    barcode_dict = generate_seq_dict(whitelist, n=1)
    linker_dict = generate_seq_dict(linker, n=2)

    fq1_list = args.fq1.split(",")
    fq2_list = args.fq2.split(",")
    # merge multiple fastq files
    if len(fq1_list) > 1:
        barcode.logger.info("merge fastq with same sample name...")
        fastq_dir = args.outdir + "/../merge_fastq"
        if not os.path.exists(fastq_dir):
            os.system('mkdir -p %s' % fastq_dir)
        fastq1_file = f"{fastq_dir}/{args.sample}_1.fq.gz"
        fastq2_file = f"{fastq_dir}/{args.sample}_2.fq.gz"
        fq1_files = " ".join(fq1_list)
        fq2_files = " ".join(fq2_list)
        fq1_cmd = f"cat {fq1_files} > {fastq1_file}"
        fq2_cmd = f"cat {fq2_files} > {fastq2_file}"
        barcode.logger.info(fq1_cmd)
        os.system(fq1_cmd)
        barcode.logger.info(fq2_cmd)
        os.system(fq2_cmd)
        barcode.logger.info("merge fastq done.")
    else:
        fastq1_file = args.fq1
        fastq2_file = args.fq2

    fh1 = xopen(fastq1_file)
    fh2 = xopen(fastq2_file)
    out_fq2 = args.outdir + '/' + args.sample + '_2.fq.gz'
    fh3 = xopen(out_fq2, 'w')

    (total_num, clean_num, no_polyT_num, lowQual_num,
     no_linker_num, no_barcode_num) = (0, 0, 0, 0, 0, 0)
    Barcode_dict = defaultdict(int)

    if args.nopolyT:
        fh1_without_polyT = xopen(args.outdir + '/noPolyT_1.fq', 'w')
        fh2_without_polyT = xopen(args.outdir + '/noPolyT_2.fq', 'w')

    if args.noLinker:
        fh1_without_linker = xopen(args.outdir + '/noLinker_1.fq', 'w')
        fh2_without_linker = xopen(args.outdir + '/noLinker_2.fq', 'w')

    bool_probe = False
    if args.probe_file and args.probe_file != 'None':
        bool_probe = True
        count_dic = genDict(dim=3)
        valid_count_dic = genDict(dim=2)
        probe_dic = read_fasta(args.probe_file)
        reads_without_probe = 0

    g1 = read_fastq(fh1)
    g2 = read_fastq(fh2)

    while True:
        try:
            (header1, seq1, qual1) = next(g1)
            (header2, seq2, qual2) = next(g2)
        except BaseException:
            break
        if total_num > 0 and total_num % 1000000 == 0:
            barcode.logger.info(
                f'processed reads: {format_number(total_num)}.'
                f'valid reads: {format_number(clean_num)}.'
            )

        total_num += 1

        # polyT filter
        if bool_T:
            polyT = seq_ranges(seq1, pattern_dict['T'])
            if no_polyT(polyT):
                no_polyT_num += 1
                if args.nopolyT:
                    fh1_without_polyT.write(
                        '@%s\n%s\n+\n%s\n' % (header1, seq1, qual1))
                    fh2_without_polyT.write(
                        '@%s\n%s\n+\n%s\n' % (header2, seq2, qual2))
                continue

        # lowQual filter
        C_U_quals_ascii = seq_ranges(
            qual1, pattern_dict['C'] + pattern_dict['U'])
        # C_U_quals_ord = [ord(q) - 33 for q in C_U_quals_ascii]
        if low_qual(C_U_quals_ascii, args.lowQual, args.lowNum):
            lowQual_num += 1
            continue

        # linker filter
        barcode_arr = [seq_ranges(seq1, [i]) for i in pattern_dict['C']]
        raw_cb = ''.join(barcode_arr)
        if bool_L:
            linker = seq_ranges(seq1, pattern_dict['L'])
            if (no_linker(linker, linker_dict)):
                no_linker_num += 1

                if args.noLinker:
                    fh1_without_linker.write(
                        '@%s\n%s\n+\n%s\n' % (header1, seq1, qual1))
                    fh2_without_linker.write(
                        '@%s\n%s\n+\n%s\n' % (header2, seq2, qual2))
                continue

        # barcode filter
        # barcode_arr = [seq_ranges(seq1, [i]) for i in pattern_dict['C']]
        # raw_cb = ''.join(barcode_arr)
        res = no_barcode(barcode_arr, barcode_dict)
        if res is True:
            no_barcode_num += 1
            continue
        elif res == "correct":
            cb = raw_cb
        else:
            cb = res

        umi = seq_ranges(seq1, pattern_dict['U'])
        Barcode_dict[cb] += 1
        clean_num += 1
        read_name_probe = 'None'

        if bool_probe:
            # valid count
            valid_count_dic[cb][umi] += 1

            # output probe UMi and read count
            find_probe = False
            for probe_name in probe_dic:
                probe_seq = probe_dic[probe_name]
                probe_seq = probe_seq.upper()
                if seq1.find(probe_seq) != -1:
                    count_dic[probe_name][cb][umi] += 1
                    read_name_probe = probe_name
                    find_probe = True
                    break

            if not find_probe:
                reads_without_probe += 1

        barcode_qual_Counter.update(C_U_quals_ascii[:C_len])
        umi_qual_Counter.update(C_U_quals_ascii[C_len:])
        C_U_base_Counter.update(raw_cb + umi)

        # new readID: @barcode_umi_old readID
        fh3.write(f'@{cb}_{umi}_{read_name_probe}_{total_num}\n{seq2}\n+\n{qual2}\n')

    fh3.close()

    # logging
    if total_num % 1000000 != 0:
        barcode.logger.info(
            f'processed reads: {format_number(total_num)}. '
            f'valid reads: {format_number(clean_num)}. '
        )

    barcode.logger.info('no_linker:',no_linker_num)
    barcode.logger.info('no_barcode:',no_barcode_num)
    barcode.logger.info('no_polyT:',no_polyT_num)

    if clean_num == 0:
        raise Exception(
            'no valid reads found! please check the --chemistry parameter.')

    if bool_probe:
        # total probe summary
        total_umi = 0
        total_valid_read = 0
        for cb in valid_count_dic:
            total_umi += len(valid_count_dic[cb])
            total_valid_read += sum(valid_count_dic[cb].values())
        barcode.logger.info("total umi:"+str(total_umi))
        barcode.logger.info("total valid read:"+str(total_valid_read))
        barcode.logger.info("reads without probe:"+str(reads_without_probe))

        # probe summary
        count_list = []
        for probe_name in probe_dic:
            UMI_count = 0
            read_count = 0
            if probe_name in count_dic:
                for cb in count_dic[probe_name]:
                    UMI_count += len(count_dic[probe_name][cb])
                    read_count += sum(count_dic[probe_name][cb].values())
            count_list.append(
                {"probe_name": probe_name, "UMI_count": UMI_count, "read_count": read_count})

        df_count = pd.DataFrame(count_list, columns=[
                                "probe_name", "read_count", "UMI_count"])

        def format_percent(x):
            x = str(round(x*100, 2))+"%"
            return x
        df_count["read_fraction"] = (
            df_count["read_count"]/total_valid_read).apply(format_percent)
        df_count["UMI_fraction"] = (
            df_count["UMI_count"]/total_umi).apply(format_percent)
        df_count.sort_values(by="UMI_count", inplace=True, ascending=False)
        df_count_file = args.outdir + '/' + args.sample + '_probe_count.tsv'
        df_count.to_csv(df_count_file, sep="\t", index=False)

    # stat
    BarcodesQ30 = sum([barcode_qual_Counter[k] for k in barcode_qual_Counter if k >= ord2chr(
        30)]) / float(sum(barcode_qual_Counter.values())) * 100
    UMIsQ30 = sum([umi_qual_Counter[k] for k in umi_qual_Counter if k >= ord2chr(
        30)]) / float(sum(umi_qual_Counter.values())) * 100

    global stat_info
    def cal_percent(x): return "{:.2%}".format((x + 0.0) / total_num)
    with open(args.outdir + '/stat.txt', 'w') as fh:
        """
        Raw Reads: %s
        Valid Reads: %s(%s)
        Q30 of Barcodes: %.2f%%
        Q30 of UMIs: %.2f%%
        """
        stat_info = stat_info % (format_number(total_num), format_number(clean_num),
                                 cal_percent(clean_num), BarcodesQ30,
                                 UMIsQ30)
        stat_info = re.sub(r'^\s+', r'', stat_info, flags=re.M)
        fh.write(stat_info)

    barcode.logger.info('fastqc ...!')
    cmd = ['fastqc', '-t', str(args.thread), '-o', args.outdir, out_fq2]
    barcode.logger.info('%s' % (' '.join(cmd)))
    subprocess.check_call(cmd)
    barcode.logger.info('fastqc done!')

    t = reporter(name='barcode', assay=args.assay, sample=args.sample,
                 stat_file=args.outdir + '/stat.txt', outdir=args.outdir + '/..')
    t.get_report()


def get_opts_barcode(parser, sub_program):
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--assay', help='assay', required=True)
    parser.add_argument('--sample', help='sample name', required=True)
    parser.add_argument('--fq1', help='read1 fq file', required=True)
    parser.add_argument('--fq2', help='read2 fq file', required=True)
    parser.add_argument(
        '--chemistry', choices=__PATTERN_DICT__.keys(), help='chemistry version')
    parser.add_argument('--pattern', help='')
    parser.add_argument('--whitelist', help='')
    parser.add_argument('--linker', help='')
    parser.add_argument('--lowQual', type=int,
                        help='max phred of base as lowQual, default=0', default=0)
    parser.add_argument(
        '--lowNum', type=int, help='max number with lowQual allowed, default=2', default=2)
    parser.add_argument('--nopolyT', action='store_true',
                        help='output nopolyT fq')
    parser.add_argument('--noLinker', action='store_true',
                        help='output noLinker fq')
    parser.add_argument('--thread', help='number of threads', default=1)
    parser.add_argument('--probe_file', help="probe fasta file")
    return parser
