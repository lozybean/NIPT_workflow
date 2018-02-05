#!/usr/bin/env python
# -*- coding: utf-8 -*- \#
"""
@author = 'liangzb'
@date = '2018/1/31 0031'

bin gc correction, cite: Fan_and_Quake
and lowess correction
"""
import argparse
import re
from functools import lru_cache

import matplotlib as mpl
import numpy as np
from NIPT_workflow.utils.bam_reader import read_bin_counts
from NIPT_workflow.utils.gc_correction import (
    read_gc_percent_file,
    gc_correct_bin_counts_expand,
    gc_correct_lowess,
)

mpl.use('Agg')
import matplotlib.pyplot as plt

re_variable = re.compile('variableStep chrom=(\w+) span=(\d+)')

gc20kBase = '/bio01/database/genome/hg19/primary23/hg19.gc20kBase.txt'


class Work(object):
    def __init__(self):
        self.gc2bin, self.bin2gc = read_gc_percent_file(self.args.gc_percent, bin_size=self.args.bin_size)
        self.bin_counts, self.gc2bin = read_bin_counts(self.args.bam_in, gc2bin=self.gc2bin,
                                                       bin_size=self.args.bin_size, min_mq=self.args.min_mq,
                                                       over_abundant=self.args.over_abundant)
        # self.gc2bin, self.bin_counts = self.filter_bins()

    @property
    @lru_cache(1)
    def args(self):
        parser = argparse.ArgumentParser(description=__doc__,
                                         formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument('-i', '--bam_in', dest='bam_in',
                            metavar='FILE', type=str, required=True,
                            help="input bam file")
        parser.add_argument('-g', '--gc_percent', dest='gc_percent',
                            metavar='FILE', type=str, default=gc20kBase,
                            help="produced by hgGcPercent,"
                                 "eg. hgGcPercent -wigOut -doGaps -win=20000 -file={gc_percent} hg19 hg19.2bit,"
                                 f"[default is {gc20kBase}]")
        parser.add_argument('-b', '--bin_size', dest='bin_size',
                            metavar='INT', type=int, default=20000,
                            help="bin size, [default is 20000]")
        parser.add_argument('-m', '--min_mq', dest='min_mq',
                            metavar='INT', type=int, default=30,
                            help="minimum mapping quality")
        parser.add_argument('-o', '--over_abundant', dest='over_abundant',
                            metavar='INT', type=int, default=3,
                            help="how many times of average depth will be recognized as "
                                 "\"OVER ABUNDANT\", [default is 3]")
        parser.add_argument('--plot_bin_count', dest='plot_bin_count',
                            metavar='FILE', type=str, default=None,
                            help="plot bin count before and after gc correction")
        parser.add_argument('--plot_bin_box', dest='plot_bin_box',
                            metavar='FILE', type=str, default=None,
                            help="box-plot bin count before and after gc correction")
        parser.add_argument('--plot_distribution', dest='plot_distribution',
                            metavar='FILE', type=str, default=None,
                            help="plot reads distribution among bins")
        args = parser.parse_args()
        return args

    @staticmethod
    def plot_bin_count_scatter(depth_dict_before, depth_dict_after, weight_dict,
                               out_file=None):
        depth_list_before = []
        depth_list_after = []
        weight_list = []
        gc_list = sorted(depth_dict_before.keys())
        for gc in gc_list:
            depth_list_before.append(depth_dict_before[gc])
            depth_list_after.append(depth_dict_after[gc])
            weight_list.append(weight_dict[gc])

        plt.figure(figsize=(10, 30))

        plt.subplot(3, 1, 1)
        plt.scatter(gc_list, depth_list_before)
        plt.xlabel('gc content(%)')
        plt.ylabel('average depth per bin')
        plt.title("before correct")

        plt.subplot(3, 1, 2)
        plt.scatter(gc_list, depth_list_after)
        plt.xlabel('gc content(%)')
        plt.ylabel('average depth per bin')
        plt.title("after correct")

        plt.subplot(3, 1, 3)
        plt.scatter(gc_list, weight_list)
        plt.xlabel('gc content(%)')
        plt.ylabel('weight')
        plt.title("weight")

        if out_file:
            plt.savefig(out_file, dpi=4000)
        else:
            plt.show()
        plt.close()

    def plot_bin_count_box(self, counts_before, counts_after, out_file=None):
        depth_list_before = []
        depth_list_after = []
        gc_list = sorted(self.gc2bin.keys())
        min_gc = min(int(gc) for gc in gc_list)
        max_gc = max(int(gc) for gc in gc_list)
        int_gc_list = np.arange(min_gc, max_gc + 0.5, 0.5)
        for gc in int_gc_list:
            list_before = []
            list_after = []
            gc_list = [i for i in self.gc2bin.keys() if gc <= i < gc + 1]
            for gc in gc_list:
                for b in self.gc2bin[gc]:
                    list_before.append(counts_before[b])
                    list_after.append(counts_after[b])
            depth_list_before.append(list_before)
            depth_list_after.append(list_after)

        plt.figure(figsize=(10, 20))

        plt.subplot(2, 1, 1)
        plt.boxplot(depth_list_before, sym='k.', positions=int_gc_list)
        plt.xlabel('gc content(%)')
        plt.ylabel('depth per bin')
        plt.title("before correct")

        plt.subplot(2, 1, 2)
        plt.boxplot(depth_list_after, sym='k.', positions=int_gc_list)
        plt.xlabel('gc content(%)')
        plt.ylabel('depth per bin')
        plt.title("after correct")

        if out_file:
            plt.savefig(out_file, dpi=4000)
        else:
            plt.show()
        plt.close()

    def plot_distribution(self, raw_counts, cor_counts, weight_dict, out_file=None):
        def plot_chrom(chrom, ind):
            dict_before = {p[1]: c for p, c in raw_counts.items() if p[0] == chrom}
            dict_after = {p[1]: c for p, c in cor_counts.items() if p[0] == chrom}
            plt.subplot(6, 4, ind)
            index = sorted(dict_before.keys())
            list_weight = []
            list_before = []
            list_after = []
            for ind in index:
                gc_content = self.bin2gc[(chrom, ind, ind + self.args.bin_size - 1)]
                weight = weight_dict[gc_content]
                list_weight.append(weight)

                list_before.append(dict_before[ind])
                list_after.append(dict_after[ind])
            epsilon = max(index)
            index = np.array(index)
            index = index * (1 / epsilon)

            plt.scatter(index, list_before, c='b', s=0.1, linewidth=0.35)
            plt.scatter(index, list_after, c='r', s=0.1, linewidth=0.35)
            plt.scatter(index, list_weight, c='y', s=0.1, linewidth=0.35)
            plt.title(chrom, size=10)
            plt.ylim(0, ymax)
            plt.xlim(0, 1)

        plt.figure(figsize=(10, 10))
        ymax = max(max(raw_counts.values()), max(cor_counts.values()))

        for i in range(1, 23):
            chrom = f'chr{i}'
            plot_chrom(chrom, i)
        plot_chrom('chrX', 23)
        plot_chrom('chrY', 24)

        plt.subplots_adjust(hspace=0.5)

        if out_file:
            plt.savefig(out_file, dpi=40)
        else:
            plt.show()
        plt.close()

    def __call__(self, *args, **kwargs):
        (gc_correct_depth, average_depth_before,
         average_depth_after, weight_dict) = gc_correct_bin_counts_expand(self.gc2bin, self.bin_counts)
        gc_correct_depth = gc_correct_lowess(self.gc2bin, gc_correct_depth)

        if self.args.plot_bin_count is not None:
            self.plot_bin_count_scatter(average_depth_before, average_depth_after, weight_dict,
                                        out_file=self.args.plot_bin_count)

        if self.args.plot_distribution is not None:
            self.plot_distribution(self.bin_counts, gc_correct_depth, weight_dict,
                                   out_file=self.args.plot_distribution)

        if self.args.plot_bin_box is not None:
            self.plot_bin_count_box(self.bin_counts, gc_correct_depth, out_file=self.args.plot_bin_box)


if __name__ == '__main__':
    work = Work()
    work()
