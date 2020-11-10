from celescope.tools.utils import format_number, log, read_barcode_file
import pandas as pd

class Barcode_index():
    def __init__(self, barcodes):
        self.barcodes = barcodes
        cell_index = 0
        index_dict = {}
        for barcode in self.barcodes:
            cell_index += 1
            index_dict[barcode] = cell_index
        self.index_dict = index_dict
        self.df_index = pd.DataFrame.from_dict(index_dict, orient='index')
        self.df_index.columns = ['cell_index']
        self.df_index.index.name = 'barcode'

    @staticmethod
    def read_index(index_file):
        df_index = pd.read_csv(index_file, sep='\t', index_col=0, dtype=object)
        df_valid = df_index[df_index['valid'] == 'True']
        return df_valid

    def write_index(self, file_name):
        """
        write index-barcode to file
        """
        self.df_index.to_csv(file_name, sep='\t')
        
    
    