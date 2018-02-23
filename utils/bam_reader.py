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
                    over_abundant=3):
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
                  min_mq=15):
    bam_file = pysam.AlignmentFile(bam_file)
    gc_bases = 0
    total_reads = 0
    total_bases = 0
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
            gc_dict[record.reference_name] += gc_base
            base_dict[record.reference_name] += len(record.seq)

    return {
        'gc_per_chrom': {chrom: gc_dict[chrom] / base_dict[chrom] * 100 for chrom in gc_dict},
        'gc_content': gc_bases / total_bases,
        'total_bases': total_bases,
        'usable_reads': total_reads,
    }


def get_reads_in_chrom(bin_counts):
    rt_dict = defaultdict(int)
    for b in bin_counts:
        rt_dict[b[0]] += bin_counts[b]
    return rt_dict


def get_ratio_in_chrom(bin_counts):
    reads_in_chrom = get_reads_in_chrom(bin_counts)
    rt_dict = {}
    total_reads = sum(reads_in_chrom.values())
    for chrom, reads_count in reads_in_chrom.items():
        rt_dict[chrom] = reads_in_chrom[chrom] / total_reads * 100
    return rt_dict
