#!/usr/bin/env python
# -*- coding: utf-8 -*- \#
"""
@author = 'liangzb'
@date = '2018/1/30 0030'

"""

import argparse
import re


def read_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-r', '--raw_bam_dir', dest='raw_bam_dir',
                        metavar='DIR', type=str, required=True,
                        help="Home dir of the results, from Ion,"
                             "eg. /analysis/output/Home/")
    parser.add_argument('-l', '--link_bam_dir', dest='link_bam_dir',
                        metavar='DIR', type=str, required=True,
                        help="directory of target link.")
    args = parser.parse_args()
    return args


def print_link(raw_bam_dir, link_bam_dir):
    re_pooling_num = re.compile('(\d+)-HZDA-HIQ-pooling')
    for pooling in raw_bam_dir.dirs('*HZDA-HIQ*'):
        pooling_num = re_pooling_num.search(pooling.namebase).group(1)
        pooling_num = int(pooling_num)
        if pooling_num > 70:
            continue
        basecaller_result = pooling / 'basecaller_results'
        for bam_file in basecaller_result.files('*.basecaller.bam'):
            if bam_file.size > 10000000 and not bam_file.namebase.startswith('nomatch'):
                link_name = link_bam_dir / bam_file.namebase.replace('IonXpress', f'Pool_{pooling_num}')
                print(f'ln -s {bam_file} {link_name}')


if __name__ == '__main__':
    args = read_args()
    print_link(args.raw_bam_dir, args.link_bam_dir)
