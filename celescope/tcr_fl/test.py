import unittest
import os
import pandas as pd
import glob
from celescope.tcr_fl.Barcode_index import Barcode_index
from celescope.tcr_fl.split_fq import *
from celescope.tcr_fl.assemble import *
from celescope.tools.utils import read_barcode_file, glob_genomeDir


class test_snp(unittest.TestCase):
    def setUp(self):
        os.chdir('/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/tcr_fl/20201103')
        self.match_dir = '/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/rna/test1'

    @unittest.skip('pass')
    def test_Barcode_index(self):
        barcodes, _nCell = read_barcode_file(self.match_dir)
        bi = Barcode_index(barcodes)
        bi.write_index('test_bi.tsv')
    
    @unittest.skip('pass')
    def test_get_nCell_barcodes(self):
        fq = 'TCR_sub_clean_2.fq.gz'
        barcodes = get_nCell_barcodes(fq, 50)
        print(barcodes)

    @unittest.skip('pass')
    def test_split_run(self):
        fq = '/SGRNJ02/RandD4/RD2019016/20201016/GOT0928_P1_T/02.cutadapt/GOT0928_P1_T_clean_2.fq.gz'
        fq_outdir = 'fastq'
        split_run(fq, fq_outdir, nCell=50)

    def test_assemble(self):
        os.chdir('/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/tcr_fl/pipe/')
        fastq_dir = 'GOT0928_P1_T/03.split_fq/fastq'
        sample = 'GOT0928_P1_T'
        outdir = 'GOT0928_P1_T/04.assemble/'
        run_assemble(sample, outdir, fastq_dir, thread=4)


if __name__ == '__main__':
    unittest.main()