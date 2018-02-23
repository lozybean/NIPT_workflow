from collections import defaultdict

import pysam
from NIPT_workflow.utils.gc_correction import filter_bins


def reads_filter(record, min_mq):
    if record.mapping_quality < min_mq:
        return True
    if record.is_supplementary:
        return True
    if record.is_secondary:
        return True
    return False


def falling_within_bin(record, bin_size):
    chrom = record.reference_name
    start = record.reference_start // bin_size * bin_size + 1
    end = start + bin_size - 1
    return chrom, start, end


def read_bin_counts(bam_file, gc2bin,
                    bin_size=20000, min_mq=15,
                    over_abundant=2):
    """
    :param bam_file: input bam file, sorted with index
    :param gc2bin: a dict like: {gc: bins}
    :param bin_size: bin size
    :param min_mq: min mapping quality
    :param over_abundant: How many times the average reads in a bin will be considered as over abundant bin.
    :return: count of bins, filtered gc2bin dict
    """
    bam_file = pysam.AlignmentFile(bam_file)
    count_bins = defaultdict(int)
    for record in bam_file:
        if reads_filter(record, min_mq):
            continue
        bin_key = falling_within_bin(record, bin_size)
        count_bins[bin_key] += 1
    gc2bin, count_bins = filter_bins(raw_counts=count_bins, gc2bin=gc2bin,
                                     over_abundant=over_abundant)
    return count_bins, gc2bin


def stat_bam_file(bam_file, accept_bins,
                  bin_size=20000,
                  min_mq=15,
                  by_window=False,
                  window_size=1000000):
    bam_file = pysam.AlignmentFile(bam_file)
    gc_bases = 0
    total_reads = 0
    total_bases = 0
    if by_window:
        gc_dict = defaultdict(lambda: defaultdict(int))
        base_dict = defaultdict(lambda: defaultdict(int))
    else:
        gc_dict = defaultdict(int)
        base_dict = defaultdict(int)
    for record in bam_file:
        if reads_filter(record, min_mq):
            continue
        bin_key = falling_within_bin(record, bin_size)
        if bin_key in accept_bins:
            total_reads += 1
            gc_base = record.seq.count('G') + record.seq.count('C')
            gc_bases += gc_base
            total_bases += len(record.seq)
            if by_window:
                resize_bin = bin_key[1] // window_size * window_size + 1
                gc_dict[record.reference_name][resize_bin] += gc_base
                base_dict[record.reference_name][resize_bin] += len(record.seq)
            else:
                gc_dict[record.reference_name] += gc_base
                base_dict[record.reference_name] += len(record.seq)
    rt_dict = {
        'gc_content': gc_bases / total_bases,
        'total_bases': total_bases,
        'usable_reads': total_reads,
    }
    if by_window:
        gc_per_window = {
            chrom: {
                window: gc_dict[chrom][window] / base_dict[chrom][window] * 100
                for window in gc_dict[chrom]
            }
            for chrom in gc_dict
        }
        rt_dict.update(gc_per_window=gc_per_window)
    else:
        gc_per_window = {chrom: gc_dict[chrom] / base_dict[chrom] * 100 for chrom in gc_dict}
        rt_dict.update(gc_per_chrom=gc_per_window)
    return rt_dict


def get_reads_in_chrom(bin_counts):
    rt_dict = defaultdict(int)
    for b in bin_counts:
        rt_dict[b[0]] += bin_counts[b]
    return rt_dict


def get_reads_in_window(bin_counts, bin_size):
    rt_dict = defaultdict(lambda: defaultdict(int))
    for b in bin_counts:
        re_size_bin = b[1] // bin_size * bin_size + 1
        rt_dict[b[0]][re_size_bin] += bin_counts[b]
    return rt_dict


def get_ratio_in_chrom(bin_counts):
    reads_in_chrom = get_reads_in_chrom(bin_counts)
    print(reads_in_chrom)
    rt_dict = {}
    total_reads = sum(reads_in_chrom.values())
    print(total_reads)
    for chrom, reads_count in reads_in_chrom.items():
        rt_dict[chrom] = reads_in_chrom[chrom] / total_reads * 100
    return rt_dict


def get_ratio_in_window(bin_counts, bin_size, min_reads_in_bins=0):
    reads_in_bins = get_reads_in_window(bin_counts, bin_size)
    rt_dict = defaultdict(dict)
    total_reads = sum(v for v in reads_in_bins.values() if v > min_reads_in_bins)
    for chrom, bins in reads_in_bins.items():
        for b in bins:
            if reads_in_bins[chrom][b] > min_reads_in_bins:
                rt_dict[chrom][b] = reads_in_bins[chrom][b] / total_reads * 100
            else:
                rt_dict[chrom][b] = 0
    return rt_dict
