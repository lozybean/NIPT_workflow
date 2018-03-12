#!/usr/bin/env python
# -*- coding: utf-8 -*- \#
"""
@author = 'liangzb'
@date = '2018/2/7 0007'

"""
import os
import sys

import django
from path import Path, getcwdu

sys.path.append(f'{Path(__file__).parent}/../server')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "server.settings")
django.setup()
from dapt.models import Sample

GENOME = '/bio01/database/genome/hg19/primary23/hg19.fasta'


def get_mapping_command(sample, mapping_dir):
    prefix = f'{sample.name}_ID_{sample.id}'
    mapping_bam = f'{mapping_dir}/{prefix}.bam'
    tempdir = f'{mapping_dir}/{prefix}.temp'
    read_group = f"'@RG\\tID:{prefix}\\tSM:{prefix}\\tLB:lib1\\tPL:IonTorrent\\tPU:unset'"
    if Path(mapping_bam).exists():
        return None
    return (f"/work/bin/samtools-1.3.1/samtools bam2fq {sample.bam_file} | "
            f"/work/bin/speedseq/bwa mem -t 8 -R {read_group} {GENOME} - | "
            f"/work/bin/speedseq/samblaster | "
            f"/work/bin/speedseq/sambamba view -S -f bam -l 0 /dev/stdin | "
            f"/work/bin/speedseq/sambamba sort -t 8 -m 8G --tmpdir {tempdir} -o {mapping_bam} "
            f"/dev/stdin && rm -R {tempdir}")


def print_command(mapping_dir):
    for sample in Sample.objects.all():
        command = get_mapping_command(sample, mapping_dir)
        if command:
            print(command)


if __name__ == '__main__':
    mapping_dir = Path(f'{getcwdu()}/../01_mapping').abspath()

    print_command(mapping_dir)
