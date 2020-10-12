import os
import glob
import sys
import argparse
import re
from collections import defaultdict
from celescope.__init__ import __CONDA__
from celescope.capture_rna.__init__ import __STEPS__, __ASSAY__
from celescope.tools.utils import merge_report, generate_sjm
from celescope.tools.utils import parse_map_col4, multi_opts, link_data


def main():

    # init
    assay = __ASSAY__
    steps = __STEPS__
    conda = __CONDA__
    app = 'celescope'

    # parser
    parser = multi_opts(assay)
    parser.add_argument('--starMem', help='starMem', default=30)
    parser.add_argument('--genomeDir', help='genome index dir', required=True)
    parser.add_argument(
        '--gtf_type',
        help='Specify attribute type in GTF annotation, default=exon',
        default='exon')
    parser.add_argument('--thread', help='thread', default=6)
    parser.add_argument('--probe_file', help="probe fasta file")
    args = parser.parse_args()

    # read args
    outdir = args.outdir
    chemistry = args.chemistry
    pattern = args.pattern
    whitelist = args.whitelist
    linker = args.linker
    lowQual = args.lowQual
    lowNum = args.lowNum
    mod = args.mod
    rm_files = args.rm_files
    outFilterMatchNmin = args.outFilterMatchNmin

    # parse mapfile
    fq_dict, match_dict = parse_map_col4(args.mapfile, None)

    # link
    link_data(outdir, fq_dict)

    # custom args
    thread = args.thread
    genomeDir = args.genomeDir
    starMem = args.starMem
    gtf_type = args.gtf_type
    probe_file = args.probe_file

    # mk log dir
    logdir = outdir + '/log'
    os.system('mkdir -p %s' % (logdir))

    # script init
    sjm_cmd = 'log_dir %s\n' % (logdir)
    sjm_order = ''
    shell_dict = defaultdict(str)

    # outdir dict
    for sample in fq_dict:
        outdir_dic = {}
        index = 0
        for step in steps:
            step_outdir = f"{outdir}/{sample}/{index:02d}.{step}"
            outdir_dic.update({step: step_outdir})
            index += 1

        # sample
        step = "sample"
        cmd = (
            f'{app} {assay} {step} '
            f'--chemistry {chemistry} '
            f'--sample {sample} --outdir {outdir_dic[step]} --assay {assay} '
        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda)
        shell_dict[sample] += cmd + '\n'
        last_step = step

        # barcode
        arr = fq_dict[sample]
        step = "barcode"
        cmd = (
            f'{app} {assay} {step} '
            f'--fq1 {arr[0]} --fq2 {arr[1]} --chemistry {chemistry} '
            f'--pattern {pattern} --whitelist {whitelist} --linker {linker} '
            f'--sample {sample} --lowQual {lowQual} --thread {thread} '
            f'--lowNum {lowNum} --outdir {outdir_dic[step]} --assay {assay} '
            f'--probe_file {probe_file} '
        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda, m=5, x=thread)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        shell_dict[sample] += cmd + '\n'
        last_step = step

        # adapt
        step = "cutadapt"
        fq = f'{outdir_dic["barcode"]}/{sample}_2.fq.gz'
        cmd = (
            f'{app} {assay} {step} '
            f'--fq {fq} --sample {sample} --outdir '
            f'{outdir_dic[step]} --assay {assay} '
        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda, m=5, x=1)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        shell_dict[sample] += cmd + '\n'
        last_step = step

        # STAR
        step = 'STAR'
        fq = f'{outdir_dic["cutadapt"]}/{sample}_clean_2.fq.gz'
        cmd = (
            f'{app} {assay} {step} '
            f'--fq {fq} --sample {sample} '
            f'--genomeDir {genomeDir} --thread {thread} ' 
            f'--outdir {outdir_dic[step]} --assay {assay} '
            f'--outFilterMatchNmin {outFilterMatchNmin} '
        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda, m=starMem, x=thread)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        shell_dict[sample] += cmd + '\n'
        last_step = step

        # featureCounts
        step = 'featureCounts'
        input = f'{outdir_dic["STAR"]}/{sample}_Aligned.sortedByCoord.out.bam'
        cmd = (
            f'{app} {assay} {step} '
            f'--input {input} --gtf_type {gtf_type} '
            f'--sample {sample} --thread {thread} --outdir {outdir_dic[step]} '
            f'--genomeDir {genomeDir} --assay {assay} '
        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda, m=8, x=thread)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        shell_dict[sample] += cmd + '\n'
        last_step = step

        # count
        step = 'count_capture_rna'
        bam = f'{outdir_dic["featureCounts"]}/{sample}_name_sorted.bam'
        cmd = (
            f'{app} {assay} {step} '
            f'--bam {bam} --sample {sample} --cells auto '
            f'--outdir {outdir_dic[step]} --assay {assay} '
            f'--match_dir {match_dict[sample]} '
            f'--genomeDir {genomeDir}'
        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda, m=8, x=thread)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        shell_dict[sample] += cmd + '\n'
        last_step = step

        # analysis
        step = 'analysis'
        matrix_file = f'{outdir_dic["count_capture_rna"]}/{sample}_matrix.tsv.gz'
        cmd = (
            f'{app} {assay} {step} '
            f'--matrix_file {matrix_file} --sample {sample} '
            f'--outdir {outdir_dic[step]} '
            f'--assay {assay} '
        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda, m=15, x=1)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        shell_dict[sample] += cmd + '\n'
        last_step = step

    # merged report
    if mod == 'sjm':
        step = 'merge_report'
        merge_report(
            fq_dict,
            steps,
            last_step,
            sjm_cmd,
            sjm_order,
            logdir,
            conda,
            outdir,
            rm_files,
        )
    if mod == 'shell':
        os.system('mkdir -p ./shell/')
        for sample in shell_dict:
            with open(f'./shell/{sample}.sh', 'w') as f:
                f.write(shell_dict[sample])


if __name__ == '__main__':
    main()
