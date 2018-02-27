from dapt.models import (
    Run, Sample, ClinicalInfo, GCContent,
    RawCoverage, Residuals, ZScore, AnalysisInfo,
    FitCoverage, FitCoverage2
)
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
admin.site.register(GCContent, AnalysisInfoAdmin)
admin.site.register(RawCoverage, AnalysisInfoAdmin)
admin.site.register(FitCoverage, AnalysisInfoAdmin)
admin.site.register(FitCoverage2, AnalysisInfoAdmin)
admin.site.register(Residuals)
admin.site.register(ZScore)
