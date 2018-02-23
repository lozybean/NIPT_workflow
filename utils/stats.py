#!/usr/bin/env python
# -*- coding: utf-8 -*- \#
"""
@author = 'liangzb'
@date = '2018/2/5 0005'

"""

from collections import defaultdict

import numpy as np
import pandas as pd
from path import Path


def get_base_name(file_name):
    return file_name.namebase.split('.')[0]


def get_stats(file_name):
    rt_dict = {}
    with open(file_name) as fp:
        for line in fp:
            if not line.startswith('#'):
                break
            key, value = line.lstrip('#').strip().split('\t')
            rt_dict[key] = float(value)
    return rt_dict


def get_one_ratio(ratio_file):
    reads_ratio = pd.read_csv(ratio_file, sep='\t', comment='#', header=0, index_col=0).T
    reads_ratio = {chrom: reads_ratio[chrom]['reads_ratio'] for chrom in reads_ratio}
    return reads_ratio


def get_one_gc(ratio_file):
    gc_content = pd.read_csv(ratio_file, sep='\t', comment='#', header=0, index_col=0).T
    gc_content = {chrom: gc_content[chrom]['gc_content'] for chrom in gc_content}
    return gc_content


def get_ratio_dict(ratio_dir, min_usable_reads=3000000, aneuploid_samples=None):
    samples = []
    gc_contents = []
    ratio_dict = defaultdict(list)
    for file in Path(ratio_dir).files():
        sample_name = get_base_name(file)
        if aneuploid_samples is not None and sample_name in aneuploid_samples:
            continue
        stats = get_stats(file)
        if stats['usable_reads'] < min_usable_reads:
            continue
        samples.append(sample_name)
        gc_contents.append(stats['gc_content'])
        reads_ratio = get_one_ratio(file)
        for chrom in reads_ratio:
            ratio_dict[chrom].append(reads_ratio[chrom])
    for chrom in ratio_dict:
        ratio_dict[chrom] = np.array(ratio_dict[chrom])
    return ratio_dict, gc_contents, samples


def ratio_fitting(ratio_dict, gc_contents):
    coefficient_dict = {}
    fit_ratio = {}
    for chrom in ratio_dict:
        coefficient = np.polyfit(gc_contents, ratio_dict[chrom], 1)
        coefficient_dict[chrom] = coefficient_dict
        fit_ratio[chrom] = np.array([coefficient[0] * gc_contents[i] + coefficient[1]
                                     for i in range(len(gc_contents))])
    return fit_ratio


def get_epsilon(ratio_dict, fit_dict):
    epsilon = {}
    for chrom in ratio_dict:
        epsilon[chrom] = ratio_dict[chrom] - fit_dict[chrom]
    return epsilon


def get_z_score(e, epsilon):
    return (e - np.mean(epsilon)) / np.std(epsilon)
