import unittest
import os
import pandas as pd
from celescope.tools.STAR import Step_mapping
from celescope.tools.utils import read_barcode_file
from .Chemistry import Chemistry

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

    @unittest.skip('pass')
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

    def test_chemistry(self):
        fq1s = [
            '/SGRNJ/DATA_PROJ/2004/20201102_9/MG_201026_1-R2010283_L3_1.fq.gz',
            '/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/test_data/smk/R2009217_L4_1.fq.gz',
            '/SGRNJ/DATA_PROJ/2003/20200727_6/R2007185_L3_1.fq.gz',
            '/SGRNJ/DATA_PROJ/2004/20201029/MG_201023_1-R2010226_L3_1.fq.gz',
        ]
        results = []
        for fq1 in fq1s:
            ch = Chemistry(fq1)
            results.append(ch.get_chemistry())
        print(results)
        assert results == ['scopeV2.1.1', 'scopeV2.1.1', 'scopeV2.0.1', 'scopeV2.2.1']            


if __name__ == '__main__':
    unittest.main()