#!/usr/bin/env python
# -*- coding: utf-8 -*- \#
"""
@author = 'liangzb'
@date = '2018/2/5 0005'

"""
import re
from collections import defaultdict

import numpy as np
from Bio.Statistics.lowess import lowess

re_variable = re.compile('variableStep chrom=(\w+) span=(\d+)')


def get_average_depth(count_dict):
    return np.mean(list(count_dict.values()))


def average_depth_in_gc(gc2bin, bin_counts, gc_content):
    bin_list = gc2bin[gc_content]
    total_depth = sum(bin_counts[b] for b in bin_list)
    return total_depth / len(bin_list)


def gc_correct_bin_counts_expand(gc2bin, raw_counts):
    rt_dict = defaultdict(int)
    average_depth_before_correct = {}
    weight_dict = {}
    average_depth_after_correct = {}
    global_average_depth = get_average_depth(raw_counts)
    for gc_content, bin_list in gc2bin.items():
        average_depth = average_depth_in_gc(gc2bin, raw_counts, gc_content)
        average_depth_before_correct[gc_content] = average_depth
        weight = global_average_depth / average_depth
        weight_dict[gc_content] = weight
        for b in bin_list:
            rt_dict[b] = raw_counts[b] * weight
        average_depth_after_correct[gc_content] = average_depth_in_gc(gc2bin, rt_dict, gc_content)
    return rt_dict, average_depth_before_correct, average_depth_after_correct, weight_dict


def gc_correct_bin_counts(gc2bin, raw_counts):
    rt_dict = defaultdict(int)
    global_average_depth = get_average_depth(raw_counts)
    for gc_content, bin_list in gc2bin.items():
        average_depth = average_depth_in_gc(gc2bin, raw_counts, gc_content)
        weight = global_average_depth / average_depth
        for b in bin_list:
            rt_dict[b] = raw_counts[b] * weight
    return rt_dict


def gc_correct_lowess(gc2bin, raw_counts):
    rt_dict = {}
    for gc_content, bin_list in gc2bin.items():
        value_list = np.array([raw_counts[b] for b in bin_list])
        if len(value_list) < 3:
            cor_value_list = value_list
        else:
            average_depth = average_depth_in_gc(gc2bin, raw_counts, gc_content)
            # key_list, value_list = zip(*sorted(bin_counts.items()))
            x = np.array(range(len(bin_list)))
            ur_loess = lowess(x, value_list)
            cor_value_list = value_list - (ur_loess - average_depth)
        for b, v in zip(bin_list, cor_value_list):
            rt_dict[b] = v
    return rt_dict


def read_gc_percent_file(file_name, bin_size):
    """
    :param file_name: get gc content file with flowing commands:
                      1. faToTwoBit hg19.fa hg19.2bit
                      2. hgGcPercent -wigOut -doGaps -win=20000 -file=gc_content_file hg19 /path/to/hg19.2bit
    :param bin_size: bin size
    :return:
    """
    gc_bin = defaultdict(list)
    bin_gc = {}
    with open(file_name) as fp:
        chrom, span = None, None
        for line in fp:
            if line.startswith('variableStep'):
                chrom, span = re_variable.search(line).groups()
                span = int(span)
                if span != bin_size:
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


def filter_bins(gc2bin, raw_counts, over_abundant):
    global_average_depth = get_average_depth(raw_counts)
    rt_gc2bin = {}
    rt_counts = {}
    for gc, bins in gc2bin.items():
        clean_bins = []
        for b in bins:
            if not raw_counts[b]:
                # filter no reads bins
                continue
            if raw_counts[b] > global_average_depth * over_abundant:
                # filter over abundant bins:
                continue
            rt_counts[b] = raw_counts[b]
            clean_bins.append(b)
        if clean_bins:
            rt_gc2bin[gc] = clean_bins
    return rt_gc2bin, rt_counts
