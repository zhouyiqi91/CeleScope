import unittest
import os
import pandas as pd
import glob
from celescope.snp.snpCalling import convert, call_all_snp, call_snp, split_bam, summary, read_index
from celescope.tools.utils import read_barcode_file, glob_genomeDir
from .analysis_snp import analysis_variant


class test_snp(unittest.TestCase):
    def setUp(self):
        os.chdir('/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/snp')
        self.genomeDir = '/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92'
        self.gene_list_file = './gene_list.tsv'
        self.index_file = './S20070818_TS/05.snpCalling/S20070818_TS_cell_index.tsv'
        self.thread = 20
        self.outdir = './S20070818_TS/05.snpCalling/'
        self.analysis_outdir = './S20070818_TS/06.analysis_snp'
        _refFlat, self.gtf, self.fasta = glob_genomeDir(self.genomeDir, fa=True)
        self.match_dir = '/SGRNJ02/RandD4/RD20051303_Panel/20200717/S20070818_ZL/'
        self.sample = 'S20070818_TS'
        self.count_file = './S20070818_TS/05.snpCalling/S20070818_TS_count.tsv'
        self.vcf_file = './S20070818_TS/05.snpCalling/S20070818_TS_anno.vcf'
        self.assay = 'snp'
        self.annovar_config = '/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/soft/annovar/annovar.config'
        self.step_analysis_variant = analysis_variant(
            self.analysis_outdir,
            self.sample,
            self.match_dir,
            self.vcf_file,
            self.index_file,
            self.assay,
            self.annovar_config,
        )

    @unittest.skip('pass')
    def test_convert(self):
        _refFlat, gtf = glob_genomeDir(self.genomeDir)
        gene_list = convert(self.gene_list_file, gtf)
        print(gene_list)

    @unittest.skip('pass')
    def test_call_all_snp(self):
        call_all_snp(self.index_file, self.outdir, self.thread, self.fasta)

    @unittest.skip('pass')
    def test_call_snp(self):
        index = 526
        call_snp(index, self.outdir, self.fasta)

    @unittest.skip('pass')
    def test_split_bam(self):
        bam = './S20070818_TS/04.featureCounts/S20070818_TS_name_sorted.bam'
        barcodes, _nCell = read_barcode_file(self.match_dir)
        gene_id_name_dic = convert(self.gene_list_file, self.gtf)
        min_query_length = 35
        split_bam(
            bam, barcodes, self.outdir,
            self.sample, gene_id_name_dic, min_query_length)
    
    @unittest.skip('pass')
    def test_summary(self):
        summary(self.index_file, self.count_file, self.outdir, self.sample)

    @unittest.skip('pass')
    def test_index(self):
        index_file = '/SGRNJ02/RandD4/RD20051303_Panel/20200929/S20070817_TS/05.snpCalling/S20070817_TS_cell_index.tsv'
        df_index, df_valid = read_index(index_file)
        index_arg = df_valid.index
        print(index_arg)
        print(list(index_arg))
    
    @unittest.skip('pass')
    def test_annovar(self):
        self.step_analysis_variant.annovar()
    
    def test_parse_annovar(self):
        from celescope.tools.utils import parse_annovar
        self.annovar_file = glob.glob(f'{self.analysis_outdir}/{self.sample}*_multianno.txt')[0]
        parse_annovar(self.annovar_file)



if __name__ == '__main__':
    unittest.main()
