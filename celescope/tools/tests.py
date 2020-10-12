import unittest
import os
import pandas as pd
from celescope.tools.STAR import Step_mapping
from celescope.tools.utils import read_barcode_file

class Tests(unittest.TestCase):
    def setUp(self):
        os.chdir('/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/rna')
        self.sample = 'test1'
        self.fq = './test1/02.cutadapt/test1_clean_2.fq.gz'
        self.genomeDir = '/SGRNJ/Public/Database/genome/homo_mus'
        self.thread = 4
        self.assay = 'rna'
        self.out_unmapped = False
        self.debug = True


    #@unittest.skip('pass')
    def test_ribo(self):
        self.outdir = f'{self.sample}/03.STAR'
        mapping = Step_mapping(
            self.sample, 
            self.outdir, 
            self.assay, 
            self.thread,
            self.fq, 
            self.genomeDir, 
            self.out_unmapped, 
            self.debug)
        mapping.STAR_map_log = f'{self.outdir}/{self.sample}_Log.final.out'  
        mapping.STAR_map_log = f'{self.outdir}/{self.sample}_Log.final.out'
        mapping.picard_region_log = f'{self.outdir}/{self.sample}_region.log'
        mapping.ribo()
        mapping.format_stat()
        mapping.report()


if __name__ == '__main__':
    unittest.main()