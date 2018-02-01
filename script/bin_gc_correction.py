#!/usr/bin/env python
# -*- coding: utf-8 -*- \#
"""
@author = 'liangzb'
@date = '2018/1/31 0031'

bin gc correction, cite: Fan_and_Quake
"""
import argparse
import re
from collections import defaultdict
from functools import lru_cache

import matplotlib as mpl
import numpy as np
import pysam

mpl.use('Agg')
import matplotlib.pyplot as plt

re_variable = re.compile('variableStep chrom=(\w+) span=(\d+)')

gc20kBase = '/mnt/analysis/NIPT_workflow/workspace/work_shell/hg19.gc20kBase.txt'


class Work(object):
    def __init__(self):
        self.gc2bin, self.bin2gc = self.read_gc_percent_file()
        self.bin_counts = self.read_bin_counts()
        self.gc2bin, self.bin_counts = self.filter_bins()

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
                            metavar='INT', type=int, default=5,
                            help="how many times of average depth will be recognized as "
                                 "\"OVER ABUNDANT\", [default is 5]")
        parser.add_argument('--plot_bin_count', dest='plot_bin_count',
                            metavar='FILE', type=str, default=None,
                            help="plot bin count before and after gc correction")
        parser.add_argument('--plot_distribution', dest='plot_distribution',
                            metavar='FILE', type=str, default=None,
                            help="plot reads distribution among bins")
        args = parser.parse_args()
        return args

    def read_gc_percent_file(self):
        gc_bin = defaultdict(list)
        bin_gc = {}
        with open(self.args.gc_percent) as fp:
            chrom, span = None, None
            for line in fp:
                if line.startswith('variableStep'):
                    chrom, span = re_variable.search(line).groups()
                    span = int(span)
                    if span != self.args.bin_size:
                        raise ValueError("bin_size do not match in gc_percent file!")
                    continue
                start, gc_content = line.strip().split('\t')
                gc_content = float(gc_content)
                if not gc_content:
                    continue
                start = int(start)
                end = start + span - 1
                gc_bin[gc_content].append((chrom, start, end))
                bin_gc[(chrom, start, end)] = gc_content
        return gc_bin, bin_gc

    def filter_bins(self):
        global_average_depth = self.get_global_average_depth()
        gc2bin = {}
        bin_counts = {}
        for gc, bins in self.gc2bin.items():
            clean_bins = []
            for b in bins:
                if not self.bin_counts[b]:
                    # filter no reads bins
                    continue
                if self.bin_counts[b] > global_average_depth * self.args.over_abundant:
                    # filter over abundant bins:
                    continue
                bin_counts[b] = self.bin_counts[b]
                clean_bins.append(b)
            if clean_bins:
                gc2bin[gc] = clean_bins
        return gc2bin, bin_counts

    def reads_filter(self, record):
        if record.mapping_quality < self.args.min_mq:
            return True
        if record.is_supplementary:
            return True
        if record.is_secondary:
            return True
        return False

    def falling_within_bin(self, record):
        chrom = record.reference_name
        start = record.reference_start // self.args.bin_size * self.args.bin_size + 1
        end = start + self.args.bin_size - 1
        return chrom, start, end

    def get_global_average_depth(self):
        return np.mean(list(self.bin_counts.values()))

    def read_bin_counts(self):
        bam_file = pysam.AlignmentFile(self.args.bam_in)
        count_bins = defaultdict(int)
        for record in bam_file:
            if self.reads_filter(record):
                continue
            bin_key = self.falling_within_bin(record)
            count_bins[bin_key] += 1
        return count_bins

    def average_depth_in_gc(self, bin_counts, gc_content):
        bin_list = self.gc2bin[gc_content]
        total_depth = sum(bin_counts[b] for b in bin_list)
        return total_depth / len(bin_list)

    def gc_correct_bin_counts(self):
        rt_dict = defaultdict(int)
        average_depth_before_correct = {}
        weight_dict = {}
        average_depth_after_correct = {}
        global_average_depth = self.get_global_average_depth()
        for gc_content, bin_list in self.gc2bin.items():
            average_depth = self.average_depth_in_gc(self.bin_counts, gc_content)
            average_depth_before_correct[gc_content] = average_depth
            weight = global_average_depth / average_depth
            weight_dict[gc_content] = weight
            for b in bin_list:
                rt_dict[b] = self.bin_counts[b] * weight
            average_depth_after_correct[gc_content] = self.average_depth_in_gc(rt_dict, gc_content)
        return rt_dict, average_depth_before_correct, average_depth_after_correct, weight_dict

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
            plt.savefig(out_file, dpi=4000)
        else:
            plt.show()
        plt.close()

    def __call__(self, *args, **kwargs):
        gc_correct_depth, average_depth_before, average_depth_after, weight_dict = self.gc_correct_bin_counts()

        if self.args.plot_bin_count is not None:
            self.plot_bin_count_scatter(average_depth_before, average_depth_after, weight_dict,
                                        out_file=self.args.plot_bin_count)

        if self.args.plot_distribution is not None:
            self.plot_distribution(self.bin_counts, gc_correct_depth, weight_dict,
                                   out_file=self.args.plot_distribution)


if __name__ == '__main__':
    work = Work()
    work()
