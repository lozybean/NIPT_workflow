import os
import sys

import django
from path import Path

sys.path.append(f'{Path(__file__).parent}/../server')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "server.settings")
django.setup()

from dapt.models import Sample, AnalysisInfo

aneuploid_samples = [
    "yang",
    "YANG",
    "yang1",
    "yang2",
    "yang3",
    "PTC",
    "ptc",
    "17HD100025",
    "17HD100114",
    "17HD120128",
    "18HD010022",
    "18HD010063R",
    "18HD010097",
    "18HD010141",
    "18HD010170",
]


def mark_control_sample():
    for sample in Sample.objects.all():
        if sample.name.startswith('17HD') or sample.name.startswith('18HD'):
            sample.is_control = False
        else:
            sample.is_control = True
        sample.save()


def mark_aneuploid_sample():
    for sample in Sample.objects.all():
        analysisinfo = getattr(sample, 'analysisinfo', AnalysisInfo())
        sample.analysisinfo = analysisinfo
        if sample.name in aneuploid_samples:
            sample.analysisinfo.aneuploid = True
        sample.analysisinfo.save()


if __name__ == '__main__':
    mark_control_sample()
    mark_aneuploid_sample()
