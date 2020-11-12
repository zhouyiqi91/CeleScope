from collections import defaultdict
import os
import pysam
from .utils import log, read_one_col
from .utils import seq_ranges, get_scope_bc, parse_pattern


class Chemistry():

    def __init__(self, fq1):
        self.fq1 = fq1
        self.nRead = 10000
        self.T4 = 0
        self.L57C = 0

    @log
    def get_chemistry(self):
        '''
        'scopeV2.0.1': 'C8L16C8L16C8L1U8T18'
        'scopeV2.1.1': 'C8L16C8L16C8L1U12T18'
        'scopeV2.2.1': 'C8L16C8L16C8L1U12T18' with 4 types of linkers
        '''
        # init
        linker_4_file, _whitelist = get_scope_bc('scopeV2.2.1')
        linker_4_list, _num = read_one_col(linker_4_file)
        linker_4_dict = defaultdict(int)
        pattern_dict = parse_pattern('C8L16C8L16C8L1U12T18')

        with pysam.FastxFile(self.fq1) as fh:
            for i in range(self.nRead):
                entry = fh.__next__()
                seq = entry.sequence
                L57C = seq[56]
                if L57C == 'C':
                    self.L57C += 1
                T4 = seq[65:69]
                if T4 == 'TTTT':
                    self.T4 += 1
                linker = seq_ranges(seq, pattern_dict=pattern_dict['L'])
                if linker in linker_4_list:
                    linker_4_dict[linker] += 1

        self.percent_T4 =  self.T4 / self.nRead 
        self.percent_L57C = self.L57C / self.nRead
        Chemistry.get_chemistry.logger.info(f'percent T4: {self.percent_T4}')
        Chemistry.get_chemistry.logger.info(f'percent L57C: {self.percent_L57C}')
        if self.percent_T4 > 0.5:
            self.chemistry = 'scopeV2.0.1'
        else:
            # V2.1.1 or V2.2.1 or failed
            self.valid_linker_type = 0
            for linker in linker_4_dict:
                linker_4_dict[linker] = linker_4_dict[linker] / self.nRead
                if linker_4_dict[linker] > 0.05:
                    self.valid_linker_type += 1
            self.linker_4_dict = linker_4_dict
            Chemistry.get_chemistry.logger.info(self.linker_4_dict)
            if self.valid_linker_type == 0:
                raise Exception('auto chemistry detection failed!')
            elif self.valid_linker_type == 1:
                self.chemistry = 'scopeV2.1.1'
            elif self.valid_linker_type < 4:
                self.chemistry = 'scopeV2.1.1'
                Chemistry.get_chemistry.logger.warning(
                    f'chemistry scopeV2.2.1 only has {self.valid_linker_type} linker types!')
            else:
                self.chemistry = 'scopeV2.2.1'
        Chemistry.get_chemistry.logger.info(f'chemistry: {self.chemistry}')
        return self.chemistry



