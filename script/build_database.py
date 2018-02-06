#!/usr/bin/env python
# -*- coding: utf-8 -*- \#
"""
@author = 'liangzb'
@date = '2018/2/5 0005'

"""

import argparse
from functools import lru_cache

import matplotlib as mpl
from NIPT_workflow.utils.stats import (
    get_ratio_dict, ratio_fitting, get_epsilon
)

mpl.use('Agg')
from matplotlib import pyplot as plt


class Work(object):

    def __init__(self):
        self.ratio_dict, self.gc_contents, self.samples = get_ratio_dict(self.args.ratio_dir,
                                                                         self.args.min_usable_reads,
                                                                         self.args.aneuploid_samples)
        self.fitting_dict = ratio_fitting(self.ratio_dict, self.gc_contents)
        self.epsilon = get_epsilon(self.ratio_dict, self.gc_contents)

    @property
    @lru_cache(1)
    def args(self):
        parser = argparse.ArgumentParser(description=__doc__,
                                         formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument('-i', '--ratio_dir', dest='ratio_dir',
                            metavar='FILE', type=str, required=True,
                            help="result directory from get_ratio_in_chrom.py")
        parser.add_argument('-m', '--min_usable_reads', dest='min_usable_reads',
                            metavar='INT', type=int, default=3000000,
                            help="minimum usable reads, [default is 3000000]")
        parser.add_argument('-f', '--aneuploid_samples', dest='aneuploid_samples',
                            metavar='STR', type=str, nargs='+',
                            help="aneuploid sample list")
        args = parser.parse_args()
        return args

    def plot_single_scatter(self, data_dict, chrom, ylabel):
        plt.scatter(self.gc_contents, data_dict[chrom])
        plt.title(chrom)
        plt.xlabel('gc content (%)')
        plt.ylabel(ylabel)
        plt.show()
        plt.close()

    def plot_all_scatter(self, data_dict):
        plt.figure(figsize=(10, 10))

        def plot_chrom(chrom, ind):
            plt.subplot(6, 4, ind)
            plt.scatter(self.gc_contents, data_dict[chrom], c='b', s=0.1, linewidth=0.35)
            plt.title(chrom, size=10)

        for i in range(1, 23):
            chrom = f'chr{i}'
            plot_chrom(chrom, i)

        plt.subplots_adjust(hspace=0.5)
        plt.show()
        plt.close()

    def plot_all_reads_ratio(self):
        self.plot_all_scatter(self.ratio_dict)

    def plot_all_epsilon(self):
        self.plot_all_scatter(self.epsilon)

    def plot_reads_ratio(self, chrom):
        self.plot_single_scatter(self.ratio_dict, chrom, ylabel='reads ratio among samples')

    def plot_epsilon_ratio(self, chrom):
        self.plot_single_scatter(self.epsilon, chrom, ylabel='epsilon among samples')
