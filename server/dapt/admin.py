from dapt.models import Run, Sample, ClinicalInfo, RawCoverage, Epsilon, ZScore, AnalysisInfo
from django.contrib import admin


# Register your models here.

class SampleAdmin(admin.ModelAdmin):
    search_fields = ['name']


class AnalysisInfoAdmin(admin.ModelAdmin):
    search_fields = ['sample__name']


admin.site.register(Run)
admin.site.register(Sample, SampleAdmin)
admin.site.register(AnalysisInfo, AnalysisInfoAdmin)
admin.site.register(ClinicalInfo)
admin.site.register(RawCoverage)
admin.site.register(Epsilon)
admin.site.register(ZScore)
