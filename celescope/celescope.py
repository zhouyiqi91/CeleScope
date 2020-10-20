#!/bin/env python
# coding=utf8

import argparse
from celescope.__init__ import __VERSION__, ASSAY_DICT


def main():
    parser = argparse.ArgumentParser(description='CeleScope')
    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version=__VERSION__)
    subparsers = parser.add_subparsers()

    # rna
    assay = 'rna'
    text = ASSAY_DICT[assay]
    subparsers_rna = subparsers.add_parser(assay, help=text, description=text)
    subparsers_rna_sub = subparsers_rna.add_subparsers()

    from celescope.tools.sample_info import sample_info, get_opts_sample
    parser_sample = subparsers_rna_sub.add_parser('sample')
    get_opts_sample(parser_sample, True)
    parser_sample.set_defaults(func=sample_info)

    from celescope.tools.barcode import barcode, get_opts_barcode
    parser_barcode = subparsers_rna_sub.add_parser('barcode')
    get_opts_barcode(parser_barcode, True)
    parser_barcode.set_defaults(func=barcode)

    from celescope.tools.cutadapt import cutadapt, get_opts_cutadapt
    parser_cutadapt = subparsers_rna_sub.add_parser('cutadapt')
    get_opts_cutadapt(parser_cutadapt, True)
    parser_cutadapt.set_defaults(func=cutadapt)

    from celescope.tools.STAR import STAR, get_opts_STAR
    parser_STAR = subparsers_rna_sub.add_parser('STAR')
    get_opts_STAR(parser_STAR, True)
    parser_STAR.set_defaults(func=STAR)

    from celescope.tools.featureCounts import featureCounts, get_opts_featureCounts
    parser_featureCounts = subparsers_rna_sub.add_parser('featureCounts')
    get_opts_featureCounts(parser_featureCounts, True)
    parser_featureCounts.set_defaults(func=featureCounts)

    from celescope.tools.count import count, get_opts_count
    parser_count = subparsers_rna_sub.add_parser('count')
    get_opts_count(parser_count, True)
    parser_count.set_defaults(func=count)

    from celescope.tools.analysis import analysis, get_opts_analysis
    parser_analysis = subparsers_rna_sub.add_parser('analysis')
    get_opts_analysis(parser_analysis, True)
    parser_analysis.set_defaults(func=analysis)

    from celescope.rna.run import run
    parser_run = subparsers_rna_sub.add_parser(
        'run', help='run all steps', conflict_handler='resolve')
    get_opts_sample(parser_run, False)
    get_opts_barcode(parser_run, False)
    get_opts_cutadapt(parser_run, False)
    get_opts_STAR(parser_run, False)
    get_opts_featureCounts(parser_run, False)
    get_opts_count(parser_run, False)
    get_opts_analysis(parser_run, False)
    parser_run.set_defaults(func=run)

    # rna_virus
    assay = 'rna_virus'
    text = ASSAY_DICT[assay]
    subparsers_rna_virus = subparsers.add_parser(
        assay, help=text, description=text)
    subparsers_rna_virus_sub = subparsers_rna_virus.add_subparsers()

    parser_sample = subparsers_rna_virus_sub.add_parser(
        'sample', description='sample infomation')
    get_opts_sample(parser_sample, True)
    parser_sample.set_defaults(func=sample_info)

    from celescope.tools.barcode import barcode, get_opts_barcode
    parser_barcode = subparsers_rna_virus_sub.add_parser('barcode')
    get_opts_barcode(parser_barcode, True)
    parser_barcode.set_defaults(func=barcode)

    from celescope.tools.cutadapt import cutadapt, get_opts_cutadapt
    parser_cutadapt = subparsers_rna_virus_sub.add_parser('cutadapt')
    get_opts_cutadapt(parser_cutadapt, True)
    parser_cutadapt.set_defaults(func=cutadapt)

    from celescope.tools.STAR import STAR, get_opts_STAR
    parser_STAR = subparsers_rna_virus_sub.add_parser('STAR')
    get_opts_STAR(parser_STAR, True)
    parser_STAR.set_defaults(func=STAR)

    from celescope.rna_virus.STAR_virus import STAR_virus, get_opts_STAR_virus
    parser_STAR_virus = subparsers_rna_virus_sub.add_parser('STAR_virus')
    get_opts_STAR_virus(parser_STAR_virus, True)
    parser_STAR_virus.set_defaults(func=STAR_virus)

    from celescope.tools.featureCounts import featureCounts, get_opts_featureCounts
    parser_featureCounts = subparsers_rna_virus_sub.add_parser('featureCounts')
    get_opts_featureCounts(parser_featureCounts, True)
    parser_featureCounts.set_defaults(func=featureCounts)

    from celescope.tools.count import count, get_opts_count
    parser_count = subparsers_rna_virus_sub.add_parser('count')
    get_opts_count(parser_count, True)
    parser_count.set_defaults(func=count)

    from celescope.rna_virus.count_virus import count_virus, get_opts_count_virus
    parser_count_virus = subparsers_rna_virus_sub.add_parser('count_virus')
    get_opts_count_virus(parser_count_virus, True)
    parser_count_virus.set_defaults(func=count_virus)

    from celescope.rna_virus.analysis_rna_virus import analysis_rna_virus, get_opts_analysis_rna_virus
    parser_analysis_rna_virus = subparsers_rna_virus_sub.add_parser(
        'analysis_rna_virus')
    get_opts_analysis_rna_virus(parser_analysis_rna_virus, True)
    parser_analysis_rna_virus.set_defaults(func=analysis_rna_virus)

    from celescope.rna_virus.run import run
    parser_run = subparsers_rna_virus_sub.add_parser(
        'run', help='run all steps', conflict_handler='resolve')
    get_opts_sample(parser_run, False)
    get_opts_barcode(parser_run, False)
    get_opts_cutadapt(parser_run, False)
    get_opts_STAR(parser_run, False)
    get_opts_STAR_virus(parser_run, False)
    get_opts_featureCounts(parser_run, False)
    get_opts_count_virus(parser_run, False)
    get_opts_analysis_rna_virus(parser_run, False)
    parser_run.set_defaults(func=run)

    # capture_virus
    assay = 'capture_virus'
    text = ASSAY_DICT[assay]
    subparsers_capture_virus = subparsers.add_parser(
        assay, help=text, description=text)
    subparsers_capture_virus_sub = subparsers_capture_virus.add_subparsers()

    parser_sample = subparsers_capture_virus_sub.add_parser('sample')
    get_opts_sample(parser_sample, True)
    parser_sample.set_defaults(func=sample_info)

    parser_barcode = subparsers_capture_virus_sub.add_parser('barcode')
    get_opts_barcode(parser_barcode, True)
    parser_barcode.set_defaults(func=barcode)

    parser_cutadapt = subparsers_capture_virus_sub.add_parser('cutadapt')
    get_opts_cutadapt(parser_cutadapt, True)
    parser_cutadapt.set_defaults(func=cutadapt)

    parser_STAR_virus = subparsers_capture_virus_sub.add_parser('STAR_virus')
    get_opts_STAR_virus(parser_STAR_virus, True)
    parser_STAR_virus.set_defaults(func=STAR_virus)

    from celescope.capture_virus.count_capture_virus import count_capture_virus, get_opts_count_capture_virus
    parser_count_capture_virus = subparsers_capture_virus_sub.add_parser(
        'count_capture_virus')
    get_opts_count_capture_virus(parser_count_capture_virus, True)
    parser_count_capture_virus.set_defaults(func=count_capture_virus)

    from celescope.capture_virus.analysis_capture_virus import analysis_capture_virus, get_opts_analysis_capture_virus
    parser_analysis_capture_virus = subparsers_capture_virus_sub.add_parser(
        'analysis_capture_virus')
    get_opts_analysis_capture_virus(parser_analysis_capture_virus, True)
    parser_analysis_capture_virus.set_defaults(func=analysis_capture_virus)

    from celescope.capture_virus.run import run
    parser_run = subparsers_capture_virus_sub.add_parser(
        'run', help='run all steps', conflict_handler='resolve')
    get_opts_sample(parser_run, False)
    get_opts_barcode(parser_run, False)
    get_opts_cutadapt(parser_run, False)
    get_opts_STAR_virus(parser_run, False)
    get_opts_count_capture_virus(parser_run, False)
    parser_run.set_defaults(func=run)

    # fusion
    assay = 'fusion'
    text = ASSAY_DICT[assay]
    subparsers_fusion = subparsers.add_parser(
        assay, help=text, description=text)
    subparsers_fusion_sub = subparsers_fusion.add_subparsers()

    parser_sample = subparsers_fusion_sub.add_parser('sample')
    get_opts_sample(parser_sample, True)
    parser_sample.set_defaults(func=sample_info)

    parser_barcode = subparsers_fusion_sub.add_parser('barcode')
    get_opts_barcode(parser_barcode, True)
    parser_barcode.set_defaults(func=barcode)

    parser_cutadapt = subparsers_fusion_sub.add_parser('cutadapt')
    get_opts_cutadapt(parser_cutadapt, True)
    parser_cutadapt.set_defaults(func=cutadapt)

    from celescope.fusion.STAR_fusion import STAR_fusion, get_opts_STAR_fusion
    parser_STAR_fusion = subparsers_fusion_sub.add_parser('STAR_fusion')
    get_opts_STAR_fusion(parser_STAR_fusion, True)
    parser_STAR_fusion.set_defaults(func=STAR_fusion)

    from celescope.fusion.count_fusion import count_fusion, get_opts_count_fusion
    parser_count_fusion = subparsers_fusion_sub.add_parser('count_fusion')
    get_opts_count_fusion(parser_count_fusion, True)
    parser_count_fusion.set_defaults(func=count_fusion)

    from celescope.fusion.run import run
    parser_run = subparsers_fusion_sub.add_parser(
        'run', help='run all steps', conflict_handler='resolve')
    get_opts_sample(parser_run, False)
    get_opts_barcode(parser_run, False)
    get_opts_cutadapt(parser_run, False)
    get_opts_STAR_fusion(parser_run, False)
    get_opts_count_fusion(parser_run, False)
    parser_run.set_defaults(func=run)

    # smk
    assay = 'smk'
    text = ASSAY_DICT[assay]
    subparsers_assay = subparsers.add_parser(
        assay, help=text, description=text)
    subparsers_assay_sub = subparsers_assay.add_subparsers()

    parser_tmp = subparsers_assay_sub.add_parser('sample')
    get_opts_sample(parser_tmp, True)
    parser_tmp.set_defaults(func=sample_info)

    parser_tmp = subparsers_assay_sub.add_parser('barcode')
    get_opts_barcode(parser_tmp, True)
    parser_tmp.set_defaults(func=barcode)

    parser_tmp = subparsers_assay_sub.add_parser('cutadapt')
    get_opts_cutadapt(parser_tmp, True)
    parser_tmp.set_defaults(func=cutadapt)

    from celescope.smk.mapping_smk import mapping_smk, get_opts_mapping_smk
    parser_tmp = subparsers_assay_sub.add_parser('mapping_smk')
    get_opts_mapping_smk(parser_tmp, True)
    parser_tmp.set_defaults(func=mapping_smk)

    from celescope.smk.count_smk import count_smk, get_opts_count_smk
    parser_tmp = subparsers_assay_sub.add_parser('count_smk')
    get_opts_count_smk(parser_tmp, True)
    parser_tmp.set_defaults(func=count_smk)

    from celescope.smk.analysis_smk import analysis_smk, get_opts_analysis_smk
    parser_tmp = subparsers_assay_sub.add_parser('analysis_smk')
    get_opts_analysis_smk(parser_tmp, True)
    parser_tmp.set_defaults(func=analysis_smk)

    from celescope.smk.run import run
    parser_tmp = subparsers_assay_sub.add_parser(
        'run', help='run all steps', conflict_handler='resolve')
    get_opts_sample(parser_tmp, False)
    get_opts_barcode(parser_tmp, False)
    get_opts_cutadapt(parser_tmp, False)
    get_opts_mapping_smk(parser_tmp, False)
    get_opts_count_smk(parser_tmp, False)
    get_opts_analysis_smk(parser_tmp, False)
    parser_tmp.set_defaults(func=run)

    # vdj
    assay = 'vdj'
    text = ASSAY_DICT[assay]
    subparsers_assay = subparsers.add_parser(
        assay, help=text, description=text)
    subparsers_assay_sub = subparsers_assay.add_subparsers()

    parser_tmp = subparsers_assay_sub.add_parser('sample')
    get_opts_sample(parser_tmp, True)
    parser_tmp.set_defaults(func=sample_info)

    parser_tmp = subparsers_assay_sub.add_parser('barcode')
    get_opts_barcode(parser_tmp, True)
    parser_tmp.set_defaults(func=barcode)

    parser_tmp = subparsers_assay_sub.add_parser('cutadapt')
    get_opts_cutadapt(parser_tmp, True)
    parser_tmp.set_defaults(func=cutadapt)

    from celescope.vdj.mapping_vdj import mapping_vdj, get_opts_mapping_vdj
    parser_tmp = subparsers_assay_sub.add_parser('mapping_vdj')
    get_opts_mapping_vdj(parser_tmp, True)
    parser_tmp.set_defaults(func=mapping_vdj)

    from celescope.vdj.count_vdj import count_vdj, get_opts_count_vdj
    parser_tmp = subparsers_assay_sub.add_parser('count_vdj')
    get_opts_count_vdj(parser_tmp, True)
    parser_tmp.set_defaults(func=count_vdj)

    from celescope.vdj.run import run
    parser_tmp = subparsers_assay_sub.add_parser(
        'run', help='run all steps', conflict_handler='resolve')
    get_opts_sample(parser_tmp, False)
    get_opts_barcode(parser_tmp, False)
    get_opts_cutadapt(parser_tmp, False)
    get_opts_mapping_vdj(parser_tmp, False)
    get_opts_count_vdj(parser_tmp, False)
    parser_tmp.set_defaults(func=run)

    '''
    # mut
    assay = 'mut'
    text = ASSAY_DICT[assay]
    subparsers_mut = subparsers.add_parser(
        assay, help=text, description=text)
    subparsers_mut_sub = subparsers_mut.add_subparsers()

    parser_sample = subparsers_mut_sub.add_parser('sample')
    get_opts_sample(parser_sample, True)
    parser_sample.set_defaults(func=sample_info)

    parser_barcode = subparsers_mut_sub.add_parser('barcode')
    get_opts_barcode(parser_barcode, True)
    parser_barcode.set_defaults(func=barcode)

    parser_cutadapt = subparsers_mut_sub.add_parser('cutadapt')
    get_opts_cutadapt(parser_cutadapt, True)
    parser_cutadapt.set_defaults(func=cutadapt)

    from celescope.mut.mapping_mut import mapping_mut, get_opts_mapping_mut
    parser_mapping_mut = subparsers_mut_sub.add_parser('mapping_mut')
    get_opts_mapping_mut(parser_mapping_mut, True)
    parser_mapping_mut.set_defaults(func=mapping_mut)

    from celescope.mut.count_mut import count_mut, get_opts_count_mut
    parser_count_mut = subparsers_mut_sub.add_parser('count_mut')
    get_opts_count_mut(parser_count_mut, True)
    parser_count_mut.set_defaults(func=count_mut)

    from celescope.fusion.run import run
    parser_run = subparsers_fusion_sub.add_parser(
        'run', help='run all steps', conflict_handler='resolve')
    get_opts_sample(parser_run, False)
    get_opts_barcode(parser_run, False)
    get_opts_cutadapt(parser_run, False)
    get_opts_STAR_fusion(parser_run, False)
    get_opts_count_fusion(parser_run, False)
    parser_run.set_defaults(func=run)
    '''

    # hla
    assay = 'hla'
    text = ASSAY_DICT[assay]
    subparsers_hla = subparsers.add_parser(
        assay, help=text, description=text)
    subparsers_hla_sub = subparsers_hla.add_subparsers()

    parser_sample = subparsers_hla_sub.add_parser('sample')
    get_opts_sample(parser_sample, True)
    parser_sample.set_defaults(func=sample_info)

    parser_barcode = subparsers_hla_sub.add_parser('barcode')
    get_opts_barcode(parser_barcode, True)
    parser_barcode.set_defaults(func=barcode)

    parser_cutadapt = subparsers_hla_sub.add_parser('cutadapt')
    get_opts_cutadapt(parser_cutadapt, True)
    parser_cutadapt.set_defaults(func=cutadapt)

    from celescope.hla.mapping_hla import mapping_hla, get_opts_mapping_hla
    parser_mapping_hla = subparsers_hla_sub.add_parser('mapping_hla')
    get_opts_mapping_hla(parser_mapping_hla, True)
    parser_mapping_hla.set_defaults(func=mapping_hla)

    '''
    from celescope.hla.count_hla import count_hla, get_opts_count_hla
    parser_count_hla = subparsers_hla_sub.add_parser('count_hla')
    get_opts_count_hla(parser_count_hla, True)
    parser_count_hla.set_defaults(func=count_hla)

    from celescope.hla.run import run
    parser_run = subparsers_fusion_sub.add_parser(
        'run', help='run all steps', conflict_handler='resolve')
    get_opts_sample(parser_run, False)
    get_opts_barcode(parser_run, False)
    get_opts_cutadapt(parser_run, False)
    get_opts_mapping_hla(parser_run, False)
    get_opts_count_hla(parser_run, False)
    parser_run.set_defaults(func=run)
    '''

    # capture_rna
    assay = 'capture_rna'
    text = ASSAY_DICT[assay]
    subparsers_capture_rna = subparsers.add_parser(assay, help=text, description=text)
    subparsers_capture_rna_sub = subparsers_capture_rna.add_subparsers()

    from celescope.tools.sample_info import sample_info, get_opts_sample
    parser_sample = subparsers_capture_rna_sub.add_parser('sample')
    get_opts_sample(parser_sample, True)
    parser_sample.set_defaults(func=sample_info)

    from celescope.tools.barcode import barcode, get_opts_barcode
    parser_barcode = subparsers_capture_rna_sub.add_parser('barcode')
    get_opts_barcode(parser_barcode, True)
    parser_barcode.set_defaults(func=barcode)

    from celescope.tools.cutadapt import cutadapt, get_opts_cutadapt
    parser_cutadapt = subparsers_capture_rna_sub.add_parser('cutadapt')
    get_opts_cutadapt(parser_cutadapt, True)
    parser_cutadapt.set_defaults(func=cutadapt)

    from celescope.tools.STAR import STAR, get_opts_STAR
    parser_STAR = subparsers_capture_rna_sub.add_parser('STAR')
    get_opts_STAR(parser_STAR, True)
    parser_STAR.set_defaults(func=STAR)

    from celescope.tools.featureCounts import featureCounts, get_opts_featureCounts
    parser_featureCounts = subparsers_capture_rna_sub.add_parser('featureCounts')
    get_opts_featureCounts(parser_featureCounts, True)
    parser_featureCounts.set_defaults(func=featureCounts)

    from celescope.capture_rna.count_capture_rna import count_capture_rna, get_opts_count_capture_rna
    parser_count_capture_rna = subparsers_capture_rna_sub.add_parser('count_capture_rna')
    get_opts_count_capture_rna(parser_count_capture_rna, True)
    parser_count_capture_rna.set_defaults(func=count_capture_rna)

    from celescope.tools.analysis import analysis, get_opts_analysis
    parser_analysis = subparsers_capture_rna_sub.add_parser('analysis')
    get_opts_analysis(parser_analysis, True)
    parser_analysis.set_defaults(func=analysis)

    # snp
    assay = 'snp'
    text = ASSAY_DICT[assay]
    subparsers_snp = subparsers.add_parser(assay, help=text, description=text)
    subparsers_snp = subparsers_snp.add_subparsers()

    from celescope.tools.sample_info import sample_info, get_opts_sample
    parser_sample = subparsers_snp.add_parser('sample')
    get_opts_sample(parser_sample, True)
    parser_sample.set_defaults(func=sample_info)

    from celescope.tools.barcode import barcode, get_opts_barcode
    parser_barcode = subparsers_snp.add_parser('barcode')
    get_opts_barcode(parser_barcode, True)
    parser_barcode.set_defaults(func=barcode)

    from celescope.tools.cutadapt import cutadapt, get_opts_cutadapt
    parser_cutadapt = subparsers_snp.add_parser('cutadapt')
    get_opts_cutadapt(parser_cutadapt, True)
    parser_cutadapt.set_defaults(func=cutadapt)

    from celescope.tools.STAR import STAR, get_opts_STAR
    parser_STAR = subparsers_snp.add_parser('STAR')
    get_opts_STAR(parser_STAR, True)
    parser_STAR.set_defaults(func=STAR)

    from celescope.tools.featureCounts import featureCounts, get_opts_featureCounts
    parser_featureCounts = subparsers_snp.add_parser('featureCounts')
    get_opts_featureCounts(parser_featureCounts, True)
    parser_featureCounts.set_defaults(func=featureCounts)

    from celescope.snp.snpCalling import snpCalling, get_opts_snpCalling
    parser_snpCalling = subparsers_snp.add_parser('snpCalling')
    get_opts_snpCalling(parser_snpCalling, True)
    parser_snpCalling.set_defaults(func=snpCalling)

    from celescope.snp.analysis_snp import analysis_snp, get_opts_analysis_snp
    parser_analysis_snp = subparsers_snp.add_parser('analysis_snp')
    get_opts_analysis_snp(parser_analysis_snp, True)
    parser_analysis_snp.set_defaults(func=analysis_snp)


    # tag
    assay = 'tag'
    text = ASSAY_DICT[assay]
    subparsers_assay = subparsers.add_parser(
        assay, help=text, description=text)
    subparsers_assay_sub = subparsers_assay.add_subparsers()

    parser_tmp = subparsers_assay_sub.add_parser('sample')
    get_opts_sample(parser_tmp, True)
    parser_tmp.set_defaults(func=sample_info)

    parser_tmp = subparsers_assay_sub.add_parser('barcode')
    get_opts_barcode(parser_tmp, True)
    parser_tmp.set_defaults(func=barcode)

    parser_tmp = subparsers_assay_sub.add_parser('cutadapt')
    get_opts_cutadapt(parser_tmp, True)
    parser_tmp.set_defaults(func=cutadapt)

    from celescope.tag.mapping_tag import mapping_tag, get_opts_mapping_tag
    parser_tmp = subparsers_assay_sub.add_parser('mapping_tag')
    get_opts_mapping_tag(parser_tmp, True)
    parser_tmp.set_defaults(func=mapping_tag)

    from celescope.tag.count_tag import count_tag, get_opts_count_tag
    parser_tmp = subparsers_assay_sub.add_parser('count_tag')
    get_opts_count_tag(parser_tmp, True)
    parser_tmp.set_defaults(func=count_tag)

    from celescope.tag.analysis_tag import analysis_tag, get_opts_analysis_tag
    parser_tmp = subparsers_assay_sub.add_parser('analysis_tag')
    get_opts_analysis_tag(parser_tmp, True)
    parser_tmp.set_defaults(func=analysis_tag)

    # citeseq
    assay = 'citeseq'
    text = ASSAY_DICT[assay]
    subparsers_assay = subparsers.add_parser(
        assay, help=text, description=text)
    subparsers_assay_sub = subparsers_assay.add_subparsers()

    parser_tmp = subparsers_assay_sub.add_parser('sample')
    get_opts_sample(parser_tmp, True)
    parser_tmp.set_defaults(func=sample_info)

    parser_tmp = subparsers_assay_sub.add_parser('barcode')
    get_opts_barcode(parser_tmp, True)
    parser_tmp.set_defaults(func=barcode)

    parser_tmp = subparsers_assay_sub.add_parser('cutadapt')
    get_opts_cutadapt(parser_tmp, True)
    parser_tmp.set_defaults(func=cutadapt)

    from celescope.tag.mapping_tag import mapping_tag, get_opts_mapping_tag
    parser_tmp = subparsers_assay_sub.add_parser('mapping_tag')
    get_opts_mapping_tag(parser_tmp, True)
    parser_tmp.set_defaults(func=mapping_tag)

    from celescope.citeseq.count_cite import count_cite, get_opts_count_cite
    parser_tmp = subparsers_assay_sub.add_parser('count_cite')
    get_opts_count_cite(parser_tmp, True)
    parser_tmp.set_defaults(func=count_cite)

    from celescope.citeseq.analysis_cite import analysis_cite, get_opts_analysis_cite
    parser_tmp = subparsers_assay_sub.add_parser('analysis_cite')
    get_opts_analysis_cite(parser_tmp, True)
    parser_tmp.set_defaults(func=analysis_cite)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
