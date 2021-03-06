#!/usr/bin/env python
# -*- coding: utf-8 -*- \#
"""
@author = 'liangzb'
@date = '2018/2/5 0005'

"""
import argparse
from functools import lru_cache

from NIPT_workflow.utils.bam_reader import (
    read_bin_counts, stat_bam_file, get_ratio_in_window
)
from NIPT_workflow.utils.gc_correction import (
    read_gc_percent_file, gc_correct_lowess, gc_correct_bin_counts
)

gc20kBase = '/bio01/database/genome/hg19/primary23/hg19.gc20kBase.txt'


class Work(object):
    def __init__(self):
        self.gc2bin, self.bin2gc = read_gc_percent_file(self.args.gc_percent, bin_size=self.args.bin_size)
        self.bin_counts, self.gc2bin = read_bin_counts(
            self.args.bam_in, gc2bin=self.gc2bin,
            bin_size=self.args.bin_size,
            min_mq=self.args.min_mq,
            over_abundant=self.args.over_abundant,
        )
        gc_correct_depth = gc_correct_bin_counts(self.gc2bin, self.bin_counts)
        self.bin_counts = gc_correct_lowess(self.gc2bin, gc_correct_depth)
        self.stats = stat_bam_file(self.args.bam_in, self.bin_counts,
                                   bin_size=self.args.bin_size,
                                   min_mq=self.args.min_mq,
                                   by_window=True,
                                   window_size=self.args.window_size)

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
        parser.add_argument('-w', '--window_size', dest='window_size',
                            metavar='INT', type=int, default=1000000,
                            help="window size, [default is 1000000]")
        parser.add_argument('-m', '--min_reads_per_window', dest='min_reads_per_window',
                            metavar='INT', type=int, default=0,
                            help="minimum reads per window")
        parser.add_argument('--min_mq', dest='min_mq',
                            metavar='INT', type=int, default=30,
                            help="minimum mapping quality")
        parser.add_argument('--over_abundant', dest='over_abundant',
                            metavar='INT', type=int, default=2,
                            help="how many times of average depth will be recognized as "
                                 "\"OVER ABUNDANT\", [default is 2]")
        parser.add_argument('-o', '--output', dest='output',
                            metavar='FILE', type=str, required=True,
                            help="output file")
        args = parser.parse_args()
        return args

    def __call__(self, *args, **kwargs):
        ratio_in_window = get_ratio_in_window(self.bin_counts, self.args.window_size)
        with open(self.args.output, 'w') as fp:
            print("#gc_content", self.stats['gc_content'] * 10000, sep='\t', file=fp)
            print('#total_bases', self.stats['total_bases'], sep='\t', file=fp)
            print('#usable_reads', self.stats['usable_reads'], sep='\t', file=fp)
            print('chrom', 'window', 'reads_ratio', 'gc_content', sep='\t', file=fp)
            for chrom in sorted(ratio_in_window.keys()):
                for window in sorted(ratio_in_window[chrom].keys()):
                    print(chrom, window, ratio_in_window[chrom][window],
                          self.stats['gc_per_window'].get(chrom, {}).get(window, 0),
                          sep='\t', file=fp)


if __name__ == '__main__':
    Work()()
