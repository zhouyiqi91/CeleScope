import os
import glob
import sys
import argparse
import re
from collections import defaultdict
from celescope.__init__ import __CONDA__
from .__init__ import __STEPS__, __ASSAY__
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
    parser.add_argument('--thread', help='thread', default=6)
    parser.add_argument("--fq_pattern", help="tag read2 pattern", required=True)
    parser.add_argument("--linker_fasta", help="linker fasta")
    parser.add_argument("--barcode_fasta", help="barcode fasta", required=True)
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
    minimum_length = args.minimum_length

    # parse mapfile
    fq_dict, match_dict = parse_map_col4(args.mapfile, None)

    # link
    link_data(outdir, fq_dict)

    # custom args
    thread = args.thread
    fq_pattern = args.fq_pattern
    linker_fasta = args.linker_fasta
    barcode_fasta = args.barcode_fasta

    # mk log dir
    logdir = outdir + '/log'
    os.system('mkdir -p %s' % (logdir))

    # script init
    sjm_cmd = 'log_dir %s\n' % (logdir)
    sjm_order = ''
    shell_dict = defaultdict(str)

    # run
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
            f'--minimum_length {minimum_length} '
        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda, m=5, x=1)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        shell_dict[sample] += cmd + '\n'
        last_step = step

        # mapping_tag
        step = 'mapping_tag'
        fq = f'{outdir_dic["cutadapt"]}/{sample}_clean_2.fq.gz'
        cmd = (
            f'{app} {assay} {step} '
            f'--sample {sample} '
            f'--outdir {outdir_dic[step]} '
            f'--assay {assay} '
            f'--fq {fq} '
            f'--fq_pattern {fq_pattern} '
            f'--barcode_fasta {barcode_fasta} '
            f'--linker_fasta {linker_fasta} '
        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda, m=5, x=1)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        shell_dict[sample] += cmd + '\n'
        last_step = step

        step = 'count_cite'
        read_count_file = f'{outdir_dic["mapping_tag"]}/{sample}_read_count.tsv'
        cmd = (
            f'{app} {assay} {step} '
            f'--sample {sample} '
            f'--outdir {outdir_dic[step]} '
            f'--assay {assay} '
            f'--match_dir {match_dict[sample]} '
            f'--read_count_file {read_count_file} '
        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda, m=5, x=1)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        shell_dict[sample] += cmd + '\n'
        last_step = step

        step = 'analysis_cite'
        citeseq_mtx = f'{outdir_dic["count_cite"]}/{sample}_citeseq.mtx.gz'
        cmd = (
            f'{app} {assay} {step} '
            f'--sample {sample} '
            f'--outdir {outdir_dic[step]} '
            f'--assay {assay} '
            f'--match_dir {match_dict[sample]} '
            f'--citeseq_mtx {citeseq_mtx} '
        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda, m=5, x=1)
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