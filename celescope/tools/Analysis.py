import os
import sys
import json
import numpy as np
import pandas as pd
import glob
from celescope.tools.report import reporter
from celescope.tools.utils import glob_genomeDir, log, parse_annovar
from celescope.tools.utils import cluster_tsne_list, marker_table, report_prepare, parse_vcf

class Analysis():

    def __init__(
        self,
        sample,
        outdir,
        assay,
        match_dir,
        step,     
    ):
        self.sample = sample
        self.outdir = outdir
        self.assay = assay
        self.match_dir = match_dir
        self.step = step

        tsne_df_file = glob.glob(f'{match_dir}/*analysis*/*tsne_coord.tsv')[0]
        marker_df_file = glob.glob(f'{match_dir}/*analysis*/*markers.tsv')[0]
        self.tsne_df = pd.read_csv(tsne_df_file, sep="\t")
        self.marker_df = pd.read_csv(marker_df_file, sep="\t")
        self.tsne_df.rename(columns={"Unnamed: 0": "barcode"}, inplace=True)
        self.cluster_tsne = cluster_tsne_list(self.tsne_df)

        if not os.path.exists(outdir):
            os.system('mkdir -p %s' % outdir)

    def add_attrs(self, *kwargs):
        for kwarg in kwargs:
            setattr(self, kwarg, kwargs[kwarg])

    @staticmethod
    def get_table(title, id, df_table):
        """
        return html code
        """
        table_dict = {}
        table_dict['title'] = title
        table_dict['table'] = df_table.to_html(
            escape=False,
            index=False,
            table_id=id,
            justify="center")
        return table_dict

    def report(self):
        t = reporter(
        name=self.step,
        assay=self.assay,
        sample=self.sample,
        outdir=self.outdir + '/..')
        t.get_report()

    @staticmethod
    def cluster_tsne_list(tsne_df, colname, show_tag=True):
        """
        tSNE_1	tSNE_2	cluster Gene_Counts
        return data list
        """
        sum_df = tsne_df.groupby([colname]).agg("count").iloc[:, 0]
        percent_df = sum_df.transform(lambda x: round(x / sum(x) * 100, 2))
        res = []
        for cluster in sorted(tsne_df[colname].unique()):
            sub_df = tsne_df[tsne_df[colname] == cluster]
            if show_tag:
                name = f"{colname} {cluster}({percent_df[cluster]}%)"
            else:
                name = f"{cluster}({percent_df[cluster]}%)"
            tSNE_1 = list(sub_df.tSNE_1)
            tSNE_2 = list(sub_df.tSNE_2)
            res.append({"name": name, "tSNE_1": tSNE_1, "tSNE_2": tSNE_2})
        return res

    def marker_table(self):
        """
        return html code
        """
        marker_df = self.marker_df.loc[:, ["cluster", "gene",
                                    "avg_logFC", "pct.1", "pct.2", "p_val_adj"]]
        marker_df["cluster"] = marker_df["cluster"].apply(lambda x: f"cluster {x}")
        marker_gene_table = marker_df.to_html(
            escape=False,
            index=False,
            table_id="marker_gene_table",
            justify="center")
        return marker_gene_table

    def report_prepare(self, **kwargs):
        json_file = self.outdir + '/../.data.json'
        if not os.path.exists(json_file):
            data = {}
        else:
            fh = open(json_file)
            data = json.load(fh)
            fh.close()

        for key in kwargs:
            data[key] = kwargs[key]

        with open(json_file, 'w') as fh:
            json.dump(data, fh)

    def run(self):
        cluster_tsne = Analysis.cluster_tsne_list(self.tsne_df, colname='cluster')
        marker_gene_table = self.marker_gene_table()
        self.report_prepare(
            cluster_tsne=cluster_tsne, 
            marker_gene_table=marker_gene_table,
        )
        self.report()