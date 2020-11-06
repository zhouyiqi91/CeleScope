from celescope.tools.report import reporter
from celescope.tools.utils import format_number, read_barcode_file, log, format_stat
import matplotlib.pyplot as plt
import os
import pandas as pd
import logging
import numpy as np
import argparse
import glob
from collections import defaultdict
import matplotlib as mpl
mpl.use('Agg')


class Count_tag():

    def __init__(
        self,
        sample,
        outdir,
        assay,
        read_count_file,
        match_dir,
        UMI_min,
        SNR_min,
        combine_cluster,
        dim,
        ):
        self.sample = sample
        self.outdir = outdir
        self.assay = assay
        self.read_count_file = read_count_file
        self.match_dir = match_dir
        self.UMI_min = UMI_min
        self.SNR_min = SNR_min
        self.combine_cluster = combine_cluster
        self.dim = int(dim)
        self.match_barcode, self.cell_total = read_barcode_file(match_dir)
        self.df_read_count = pd.read_csv(read_count_file, sep="\t", index_col=0)
        self.tsne_file = glob.glob(f'{match_dir}/*analysis/*tsne_coord.tsv')[0]
        
        if not os.path.exists(outdir):
            os.system('mkdir -p %s' % outdir)


    @staticmethod
    def get_UMI(row):
        return row.sum()

    @staticmethod
    def get_UMI_min(df_cell_UMI, UMI_min):
        if UMI_min == "auto":
            UMI_min1 = np.percentile(df_cell_UMI.sum(axis=1), 5)
            UMI_min2 = np.median(df_cell_UMI.sum(axis=1)) / 10
            UMI_min = int(min(UMI_min1, UMI_min2))
            UMI_min = max(UMI_min, 1)
            return UMI_min
        else:
            return int(UMI_min)

    @staticmethod
    def get_SNR(row, dim):
        row_sorted = sorted(row, reverse=True)
        noise = row_sorted[dim]
        signal = row_sorted[dim - 1]
        if noise == 0:
            return np.inf
        return float(signal) / noise

    @staticmethod
    def get_SNR_min(df_cell_UMI, dim, SNR_min, UMI_min):
        UMIs = df_cell_UMI.apply(Count_tag.get_UMI, axis=1)
        df_valid_cell_UMI = df_cell_UMI[UMIs >= UMI_min]
        if SNR_min == "auto":
            SNRs = df_valid_cell_UMI.apply(Count_tag.get_SNR, dim=dim, axis=1)
            if np.median(SNRs) == np.inf:
                return 10
            return max(np.median(SNRs) / 10, 2)
        else:
            return float(SNR_min)

    @staticmethod
    def tag_type(row, UMI_min, SNR_min, dim):
        SNR = Count_tag.get_SNR(row, dim)
        UMI = Count_tag.get_UMI(row)
        if UMI < UMI_min:
            return "Undetermined"
        if SNR < SNR_min:
            return "Multiplet"
        # get tag
        signal_tags = sorted(row.sort_values(ascending=False).index[0:dim])
        signal_tags_str = "_".join(signal_tags)
        return signal_tags_str


    def write_and_plot(self, df, column_name, count_file, plot_file):
        df_count = df.groupby(["tag", column_name]).size().unstack()
        df_count.fillna(0, inplace=True)
        df_count.to_csv(count_file, sep="\t")
        df_percent = df_count / df_count.sum()
        df_plot = df_percent.stack().reset_index()
        df_plot.rename({0: "percent"}, axis=1, inplace=True)

        # plot
        colors = list(mpl.colors.cnames.keys())
        fig, ax = plt.subplots(figsize=(20, 10))
        types = df_plot["tag"].drop_duplicates()
        margin_bottom = np.zeros(len(df_plot[column_name].drop_duplicates()))

        for num, type in enumerate(types):
            values = list(df_plot.loc[df_plot["tag"] == type, "percent"])
            df_plot[df_plot['tag'] == type].plot.bar(
                x=column_name, y='percent', ax=ax, stacked=True,
                bottom=margin_bottom, label=type, color=colors[num * 3 + 1])
            margin_bottom += values
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.title("tag fraction")
        fig.savefig(plot_file)

    @log
    def run(self):
        stats = pd.Series()
        outdir = self.outdir
        sample = self.sample
        UMI_tag_file = f'{outdir}/{sample}_umi_tag.tsv'
        tsne_tag_file = f'{outdir}/{sample}_tsne_tag.tsv'
        cluster_count_file = f'{outdir}/{sample}_cluster_count.tsv'
        cluster_plot = f'{outdir}/{sample}_cluster_plot.pdf'
        if self.combine_cluster:
            combine_cluster_count_file = f'{outdir}/{sample}_combine_cluster_count.tsv'
            combine_cluster_plot = f'{outdir}/{sample}_combine_cluster_plot.pdf'

        mapped_read = self.df_read_count['read_count'].sum()

        # in cell
        df_read_count_in_cell = self.df_read_count[self.df_read_count.index.isin(self.match_barcode)]
        mapped_read_in_cell = int(df_read_count_in_cell['read_count'].sum())
        stats = stats.append(pd.Series(
            format_stat(mapped_read_in_cell, mapped_read),
            index=['Mapped Reads in Cells']
        ))

        # UMI
        tag_name = df_read_count_in_cell.columns[0]
        print(tag_name)
        df_UMI_in_cell = df_read_count_in_cell.reset_index().groupby([
            'barcode', tag_name]).agg({'UMI': 'count'})
        df_UMI_in_cell = df_UMI_in_cell.reset_index()
        df_UMI_in_cell = df_UMI_in_cell.pivot(
            index='barcode', columns=tag_name, values='UMI')
        df_cell = pd.DataFrame(index=self.match_barcode)
        df_UMI_cell = pd.merge(
            df_cell,
            df_UMI_in_cell,
            how="left",
            left_index=True,
            right_index=True)

        # fillna
        df_UMI_cell.fillna(0, inplace=True)
        df_UMI_cell = df_UMI_cell.astype(int)

        # UMI
        UMIs = df_UMI_cell.apply(sum, axis=1)
        median = round(np.median(UMIs), 2)
        mean = round(np.mean(UMIs), 2)
        stats = stats.append(pd.Series(
            str(median),
            index=['Median UMI per Cell']
        ))

        stats = stats.append(pd.Series(
            str(mean),
            index=['Mean UMI per Cell']
        ))

        UMI_min = Count_tag.get_UMI_min(df_UMI_cell, self.UMI_min)
        Count_tag.run.logger.info(f'UMI_min: {UMI_min}')
        SNR_min = Count_tag.get_SNR_min(df_UMI_cell, self.dim, self.SNR_min, UMI_min)
        Count_tag.run.logger.info(f'SNR_min: {SNR_min}')
        df_UMI_cell["tag"] = df_UMI_cell.apply(
            Count_tag.tag_type, UMI_min=UMI_min, SNR_min=SNR_min, dim=self.dim, axis=1)
        df_UMI_cell.to_csv(UMI_tag_file, sep="\t")

        df_tsne = pd.read_csv(self.tsne_file, sep="\t", index_col=0)
        df_tsne_tag = pd.merge(
            df_tsne,
            df_UMI_cell,
            how="left",
            left_index=True,
            right_index=True)

        if self.combine_cluster:
            df_combine_cluster = pd.read_csv(
                self.combine_cluster, sep="\t", header=None)
            df_combine_cluster.columns = ["cluster", "combine_cluster"]
            df_tsne_combine_cluster_tag = pd.merge(
                df_tsne_tag, df_combine_cluster,
                on=["cluster"], how="left", left_index=True).set_index(df_tsne_tag.index)
            df_tsne_combine_cluster_tag.to_csv(tsne_tag_file, sep="\t")
        else:
            df_tsne_tag.to_csv(tsne_tag_file, sep="\t")

        self.write_and_plot(
            df=df_tsne_tag, column_name="cluster", count_file=cluster_count_file,
            plot_file=cluster_plot
        )

        if self.combine_cluster:
            self.write_and_plot(
                df=df_tsne_combine_cluster_tag,
                column_name="combine_cluster",
                count_file=combine_cluster_count_file,
                plot_file=combine_cluster_plot
            )

        df_tag_count = df_UMI_cell["tag"].value_counts().reset_index()
        df_tag_count.columns = ["item", "count"]
        for index, row in df_tag_count.iterrows():
            stats = stats.append(pd.Series(
                format_stat(row['count'], self.cell_total),
                index=[row['item'] + ' Cells']
        ))
        self.stats = stats

    def report(self):

        self.stat_file = f'{self.outdir}/stat.txt'
        self.stats.to_csv(self.stat_file, sep=':', header=False)
        t = reporter(
        name='count_tag', 
        assay=self.assay, 
        sample=self.sample,
        stat_file=self.stat_file, 
        outdir=self.outdir + '/..')
        t.get_report()
