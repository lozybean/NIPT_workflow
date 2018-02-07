#!/usr/bin/env python
# -*- coding: utf-8 -*- \#
"""
@author = 'liangzb'
@date = '2018/1/30 0030'

"""

import os
import sys

import django
from path import Path, getcwdu

sys.path.append(f'{Path(__file__).parent}/../server')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "server.settings")
django.setup()

from dapt.models import Sample


def print_link(link_bam_dir):
    for sample in Sample.objects.all():
        link_name = (f'{link_bam_dir}/{sample.name}_{sample.run.pool_number}_'
                     f'{sample.run.pool_lab}_{sample.seq_barcode}.bam')

        print(f'ln -s {sample.bam_file} {link_name}')


if __name__ == '__main__':
    link_bam_dir = Path(f'{getcwdu()}/../00_raw_bam').abspath()
    link_bam_dir.mkdir_p()
    print_link(link_bam_dir)
