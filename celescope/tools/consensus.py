import pysam
import numpy as np
import subprocess
from xopen import xopen
from collections import defaultdict
from itertools import groupby

from celescope.tools.Step import Step, s_common
import celescope.tools.utils as utils


@utils.add_log
def sort_fastq(fq, fq_tmp_file, outdir):
    tmp_dir = f'{outdir}/tmp'
    cmd = (
        f'mkdir {tmp_dir};'
        f'less {fq} | paste - - - - | sort -T {tmp_dir} -k1,1 -t " " | tr "\t" "\n" > {fq_tmp_file};'
    )
    subprocess.check_call(cmd, shell=True)


@utils.add_log
def sorted_dumb_consensus(fq, outfile, threshold):
    '''
    consensus read in name-sorted fastq
    output (barcode,umi) consensus fastq
    '''
    read_list = []
    n_umi = 0
    total_ambiguous_base_n = 0
    length_list = []
    out_h = xopen(outfile, 'w')

    def keyfunc(read):
        attr = read.name.split('_')
        return (attr[0], attr[1])

    with pysam.FastxFile(fq) as fh:
        for (barcode, umi), g in groupby(fh, key=keyfunc):
            read_list = []
            for read in g:
                read_list.append([read.sequence, read.quality])
            consensus_seq, consensus_qual, ambiguous_base_n, con_len = dumb_consensus(
                read_list, threshold=threshold, ambiguous="N")
            n_umi += 1
            prefix = "_".join([barcode, umi])
            read_name = f'{prefix}_{n_umi}'
            out_h.write(utils.fastq_line(read_name, consensus_seq, consensus_qual))
            if n_umi % 10000 == 0:
                sorted_dumb_consensus.logger.info(f'{n_umi} UMI done.')
            total_ambiguous_base_n += ambiguous_base_n
            length_list.append(con_len)
    
    out_h.close()
    return n_umi, total_ambiguous_base_n, length_list


@utils.add_log
def wrap_consensus(fq, outdir, sample, threshold):
    fq_tmp_file = f'{outdir}/{sample}_sorted.fq.tmp'
    sort_fastq(fq, fq_tmp_file, outdir)
    outfile = f'{outdir}/{sample}_consensus.fq'
    n, total_ambiguous_base_n, length_list = sorted_dumb_consensus(
        fq=fq_tmp_file, outfile=outfile, threshold=threshold)
    return outfile, n, total_ambiguous_base_n, length_list


def dumb_consensus(read_list, threshold=0.5, ambiguous='N', default_qual='F'):
    '''
    This is similar to biopython dumb_consensus.
    It will just go through the sequence residue by residue and count up the number of each type
    of residue (ie. A or G or T or C for DNA) in all sequences in the
    alignment. If the percentage of the most common residue type is
    greater then the passed threshold, then we will add that residue type,
    otherwise an ambiguous character will be added.
    elements of read_list: [entry.sequence,entry.quality]
    '''

    con_len = get_read_length(read_list, threshold=threshold)
    consensus_seq = ""
    consensus_qual = ""
    ambiguous_base_n = 0
    for n in range(con_len):
        atom_dict = defaultdict(int)
        quality_dict = defaultdict(int)
        num_atoms = 0
        for read in read_list:
            # make sure we haven't run past the end of any sequences
            # if they are of different lengths
            sequence = read[0]
            quality = read[1]
            if n < len(sequence):
                atom = sequence[n]
                atom_dict[atom] += 1
                num_atoms = num_atoms + 1

                base_qual = quality[n]
                quality_dict[base_qual] += 1

        consensus_atom = ambiguous
        for atom in atom_dict:
            if atom_dict[atom] >= num_atoms * threshold:
                consensus_atom = atom
                break
        if consensus_atom == ambiguous:
            ambiguous_base_n += 1
        consensus_seq += consensus_atom

        max_freq_qual = 0
        consensus_base_qual = default_qual
        for base_qual in quality_dict:
            if quality_dict[base_qual] > max_freq_qual:
                max_freq_qual = quality_dict[base_qual]
                consensus_base_qual = base_qual

        consensus_qual += consensus_base_qual
    return consensus_seq, consensus_qual, ambiguous_base_n, con_len


def get_read_length(read_list, threshold=0.5):
    '''
    compute read_length from read_list. 
    length = max length with read fraction >= threshold
    elements of read_list: [entry.sequence,entry.quality]
    '''
    
    n_read = len(read_list)
    length_dict = defaultdict(int)
    for read in read_list:
        length = len(read[0])
        length_dict[length] += 1
    for length in length_dict:
        length_dict[length] = length_dict[length] / n_read

    fraction = 0
    for length in sorted(length_dict.keys(),reverse=True):
        fraction += length_dict[length]
        if fraction >= threshold:
            return length

@utils.add_log
def consensus(args):

    step_name = "consensus"
    step = Step(args, step_name)

    if args.not_consensus:
        consensus.logger.warning("not_consensus!")
        return
    sample = args.sample
    outdir = args.outdir
    assay = args.assay
    fq = args.fq
    threshold = float(args.threshold)

    outfile, n, total_ambiguous_base_n, length_list = wrap_consensus(fq, outdir, sample, threshold)

    # metrics
    metrics = {}
    metrics["UMI Counts"] = n
    metrics["Mean UMI Length"] = np.mean(length_list)
    metrics["Ambiguous Base Counts"] = total_ambiguous_base_n    
    utils.format_metrics(metrics)

    ratios = {}
    ratios["Ambiguous Base Counts Ratio"] = total_ambiguous_base_n / sum(length_list)
    utils.format_ratios(ratios)

    # stat file
    stat_file = f'{outdir}/stat.txt'
    with open(stat_file, 'w') as stat_h:
        stat_str = (
            f'UMI Counts: {metrics["UMI Counts"]}\n'
            f'Mean UMI Length: {metrics["Mean UMI Length"]}\n'
            f'Ambiguous Base Counts: {metrics["Ambiguous Base Counts"]}({ratios["Ambiguous Base Counts Ratio"]}%)\n'
        )
        stat_h.write(stat_str)

    step.clean_up()



def get_opts_consensus(parser, sub_program):
    parser.add_argument("--threshold", help='valid base threshold', default=0.5)
    parser.add_argument("--not_consensus", help="input fastq is not consensus", action='store_true')
    if sub_program:
        s_common(parser)
        parser.add_argument("--fq", required=True)