import os
import glob
import sys
import argparse
import re
from collections import defaultdict
from celescope.__init__ import __CONDA__
from celescope.hla.__init__ import __STEPS__, __ASSAY__
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
    parser.add_argument('--genomeDir', help='genome index dir', required=True)
    parser.add_argument('--starMem', help='starMem', default=30)
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

    # parse mapfile
    fq_dict, match_dict = parse_map_col4(args.mapfile, None)

    # link
    link_data(outdir, fq_dict)

    # custom args
    thread = args.thread
    genomeDir = args.genomeDir
    starMem = args.starMem

    # mk log dir
    logdir = outdir + '/log'
    os.system('mkdir -p %s' % (logdir))

    # script init
    sjm_cmd = 'log_dir %s\n' % (logdir)
    sjm_order = ''
    shell = ''
    

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
        shell += cmd + '\n'
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
        shell += cmd + '\n'
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
        shell += cmd + '\n'
        last_step = step

        # STAR 
        
        step = "STAR"
        fq = f'{outdir_dic["cutadapt"]}/{sample}_clean_2.fq.gz'
        cmd = (
            f'{app} {assay} {step} '
            f'--fq {fq} --sample {sample} '
            f'--genomeDir {genomeDir} --thread {thread} ' 
            f'--outdir {outdir_dic[step]} --assay {assay} '
        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda, m=starMem, x=thread)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        shell += cmd + '\n'
        last_step = step
        

        
        # mapping_hla
        step = 'mapping_hla'
        bam = f'{outdir_dic["STAR"]}/{sample}_Aligned.sortedByCoord.out.bam'
        cmd = (
            f'{app} {assay} {step} '
            f'--bam {bam} --sample {sample} '
            f'--thread {thread} '
            f'--match_dir {match_dict[sample]} '
            f'--outdir {outdir_dic[step]} --assay {assay} '
        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda, m=30, x=thread)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        shell += cmd + '\n'
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
        with open(f'./shell/{sample}.sh', 'w') as f:
            f.write(shell)


if __name__ == '__main__':
    main()
