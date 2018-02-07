#!/usr/bin/env python
# -*- coding: utf-8 -*- \#
"""
@author = 'liangzb'
@date = '2018/1/24 0024'

"""

import json
import os
import re
import sys

import django
from path import Path

sys.path.append(f'{Path(__file__).parent}/../server')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "server.settings")
django.setup()

from dapt.models import Run, Sample


def insert_a_run(dirname):
    re_name = re.compile('Auto_\S+-(\d+)-HZDA-HIQ-pooling-(\d+)-(\d+)*')
    name_search = re_name.search(dirname.namebase)
    if not name_search:
        return None
    pool_number, pool_lab, seq_date = re_name.search(dirname.namebase).groups()
    home_dir = dirname
    r = Run(pool_number=pool_number, pool_lab=pool_lab, home_dir=home_dir)
    r.save()
    return r


def insert_a_sample(run, sample_name, infos):
    barcode = int(infos['barcodes'][0].split('_')[-1])
    original_bam = Path(f"{run.home_dir}/basecaller_results/{infos['barcodes'][0]}_rawlib.basecaller.bam")
    if not original_bam.exists():
        return None
    sample = Sample(name=sample_name, seq_barcode=barcode, run=run, bam_file=original_bam)
    sample.save()
    return sample


def insert_samples(run):
    ion_params = Path(run.home_dir) / 'ion_params_00.json'
    sample_list = []
    with open(ion_params) as fp:
        ion_params = json.load(fp)
    for sample_name, infos in ion_params['experimentAnalysisSettings']['barcodedSamples'].items():
        sample = insert_a_sample(run, sample_name, infos)
        sample_list.append(sample)
    return sample_list


def insert_runs(base_dir):
    base_dir = Path(base_dir)
    re_end = re.compile('.*_tn_\d+$')
    for dir_name in base_dir.dirs():
        if re_end.match(dir_name):
            continue
        run = insert_a_run(dir_name)
        if run is None:
            continue
        samples = insert_samples(run)
        if set(samples) == {None}:
            run.delete()


if __name__ == '__main__':
    insert_runs('/mnt/rawdata/NIPT_backup/results/Home/')
