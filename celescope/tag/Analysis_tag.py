import os
import sys
import json
import numpy as np
import pandas as pd
import glob
from celescope.tools.report import reporter
from celescope.tools.utils import *
from celescope.tools.Analysis import Analysis

class Analysis_tag(Analysis):

    def run(self, tsne_tag_file):
        self.cluster_tsne = self.get_cluster_tsne(colname='cluster')
        self.tsne_tag_df = pd.read_csv(tsne_tag_file, sep="\t", index_col=0)
        self.feature_tsne = self.get_cluster_tsne(colname='tag', show_tag=False, dfname='tsne_tag_df')
        self.marker_gene_table = self.marker_table()






