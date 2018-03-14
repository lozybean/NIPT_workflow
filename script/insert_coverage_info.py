import os
import re
import sys
from collections import defaultdict

import django
import numpy as np
from NIPT_workflow.utils.bam_reader import read_bin_counts, get_ratio_in_chrom, stat_bam_file
from NIPT_workflow.utils.gc_correction import read_gc_percent_file
from NIPT_workflow.utils.stats import get_one_ratio, get_stats, get_one_gc
from path import Path, getcwdu

sys.path.append(f'{Path(__file__).parent}/../server')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "server.settings")
django.setup()

from dapt.models import Sample, AnalysisInfo, RawCoverage, RawCoverage2, GCContent, FitCoverage, FitCoverage2

re_id = re.compile('ID_(\d+)')
gc20kBase = '/bio01/database/genome/hg19/primary23/hg19.gc20kBase.txt'
chrom_normal = [f'chr{i}' for i in range(1, 23)]
chrom_list = chrom_normal + ['chrX', 'chrY']


def get_coverage(sample, bam_file):
    sample.bam_file = bam_file
    sample.save()

    gc2bin, bin2gc = read_gc_percent_file(gc20kBase, bin_size=20000)
    bin_counts, gc2bin = read_bin_counts(sample.bam_file, gc2bin=gc2bin,
                                         bin_size=20000, min_mq=15, over_abundant=2)
    stats = stat_bam_file(sample.bam_file, bin_counts, bin_size=20000, min_mq=15)
    ratio_in_chrom = get_ratio_in_chrom(bin_counts)

    sample.analysis_info = AnalysisInfo(sample=sample,
                                        gc_percent=stats['gc_content'],
                                        usable_reads=stats['usable_reads'])
    sample.rawcoverage = RawCoverage(**ratio_in_chrom, sample=sample)

    sample.save()


def get_sample(ratio_file):
    sample_id = re_id.search(Path(ratio_file).namebase).group(1)
    sample = Sample.objects.get(pk=sample_id)
    return sample


def read_ratio(sample, ratio_file):
    ratio_dict = get_one_ratio(ratio_file)
    coverage = getattr(sample, 'rawcoverage', RawCoverage())
    for chrom in chrom_list:
        value = ratio_dict.get(chrom, 0)
        setattr(coverage, chrom, value)
    sample.rawcoverage = coverage
    sample.rawcoverage.save()


def read_ratio2(sample, ratio_file):
    ratio_dict = get_one_ratio(ratio_file)
    coverage = getattr(sample, 'rawcoverage2', RawCoverage2())
    for chrom in chrom_list:
        value = ratio_dict.get(chrom, 0)
        setattr(coverage, chrom, value)
    sample.rawcoverage2 = coverage
    sample.rawcoverage2.save()


def read_gc(sample, ratio_file):
    gc_dict = get_one_gc(ratio_file)
    gc_content = getattr(sample, 'gccontent', GCContent())
    for chrom in chrom_list:
        value = gc_dict.get(chrom, 0)
        setattr(gc_content, chrom, value)
    sample.gccontent = gc_content
    sample.gccontent.save()


def read_stat(sample, ratio_file):
    analysis_info = getattr(sample, 'analysisinfo', AnalysisInfo())
    stat = get_stats(ratio_file)
    analysis_info.gc_content = stat['gc_content']
    analysis_info.usable_reads = stat['usable_reads']
    sample.analysisinfo = analysis_info
    sample.analysisinfo.save()


def read_all_ratio(ratio_dir):
    for file in Path(ratio_dir).files('*.ratio.txt'):
        sample = get_sample(file)
        read_ratio(sample, file)
        read_gc(sample, file)
        read_stat(sample, file)


def read_all_ratio2(ratio_dir):
    for file in Path(ratio_dir).files('*.ratio.txt'):
        sample = get_sample(file)
        read_ratio2(sample, file)


def fitting_all_ratio():
    ratio_dict = defaultdict(list)
    gc_content_list = []
    sample_list = []
    for sample in Sample.objects.filter(is_control=False, analysisinfo__aneuploid=False,
                                        analysisinfo__usable_reads__gt=3000000,
                                        is_abnormal=False):
        sample_list.append(sample)
        gc_content_list.append(sample.analysisinfo.gc_content)
        for chrom in chrom_list:
            ratio_dict[chrom].append(getattr(sample.rawcoverage, chrom))
    coefficient_dict = {}
    for chrom in ratio_dict:
        coefficient_dict[chrom] = np.polyfit(gc_content_list, ratio_dict[chrom], 1)
    for sample in Sample.objects.all():
        fit_ratio = {}
        for chrom in chrom_list:
            sample_gc_chrom = sample.analysisinfo.gc_content
            fit_ratio[chrom] = coefficient_dict[chrom][0] * sample_gc_chrom + coefficient_dict[chrom][1]
        coverage = getattr(sample, 'fitcoverage', FitCoverage())
        for chrom in fit_ratio:
            setattr(coverage, chrom, fit_ratio[chrom])
        sample.fitcoverage = coverage
        sample.fitcoverage.save()


def fitting_all_ratio2():
    ratio_dict = defaultdict(list)
    gc_dict = defaultdict(list)
    sample_list = []
    for sample in Sample.objects.filter(is_control=False, analysisinfo__aneuploid=False,
                                        analysisinfo__usable_reads__gt=3000000,
                                        is_abnormal=False):
        sample_list.append(sample)
        for chrom in chrom_list:
            ratio_dict[chrom].append(getattr(sample.rawcoverage, chrom))
        for chrom in chrom_list:
            gc_dict[chrom].append(getattr(sample.gccontent, chrom))
    coefficient_dict = {}
    for chrom in ratio_dict:
        coefficient_dict[chrom] = np.polyfit(gc_dict[chrom], ratio_dict[chrom], 1)
    for sample in Sample.objects.all():
        fit_ratio = {}
        for chrom in chrom_list:
            sample_gc_chrom = getattr(sample.gccontent, chrom)
            fit_ratio[chrom] = coefficient_dict[chrom][0] * sample_gc_chrom + coefficient_dict[chrom][1]
        coverage = getattr(sample, 'fitcoverage2', FitCoverage2())
        for chrom in fit_ratio:
            setattr(coverage, chrom, fit_ratio[chrom])
        sample.fitcoverage2 = coverage
        sample.fitcoverage2.save()


if __name__ == '__main__':
    ratio_dir = Path(f'{getcwdu()}/../03_reads_ratio').abspath()
    read_all_ratio(ratio_dir)
    fitting_all_ratio()
    fitting_all_ratio2()
    # ratio_dir = Path(f'{getcwdu()}/../03_reads_ratio_without_gc_correct').abspath()
    # read_all_ratio2(ratio_dir)
