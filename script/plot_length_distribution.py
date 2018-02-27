#!/usr/bin/env python
# -*- coding: utf-8 -*- \#
"""
@author = 'liangzb'
@date = '2018/2/24 0024'

"""

import argparse
from functools import lru_cache

import matplotlib as mpl
from NIPT_workflow.utils.bam_reader import (
    read_length_distribution,
)

mpl.use('Agg')


class Work(object):
    def __init__(self):
        self.distribution = read_length_distribution(self.args.bam_in, min_mq=self.args.min_mq)

    @property
    @lru_cache(1)
    def args(self):
        parser = argparse.ArgumentParser(description=__doc__,
                                         formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument('-i', '--bam_in', dest='bam_in',
                            metavar='FILE', type=str, required=True,
                            help="input bam file")
        parser.add_argument('-m', '--min_mq', dest='min_mq',
                            metavar='INT', type=int, default=30,
                            help="minimum mapping quality, [default is 30]")
        parser.add_argument('-t', '--interval', dest='interval',
                            metavar='INT', type=int, default=5,
                            help="set the interval, [default is 5]")
        parser.add_argument('-o', '--output', dest='output',
                            metavar='FILE', type=str, required=True,
                            help="output file")
        args = parser.parse_args()
        return args

    def __call__(self, *args, **kwargs):
        with open(self.args.output, 'w') as fp:
            print('length', 'reads', sep='\t', file=fp)
            for length in sorted(self.distribution.keys()):
                # for key, value in self.distribution.items():
                print(length, self.distribution[length], sep='\t', file=fp)


if __name__ == '__main__':
    Work()()
