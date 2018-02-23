from dapt.models import Run, Sample, ClinicalInfo, RawCoverage, Residuals, ZScore, AnalysisInfo
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
admin.site.register(RawCoverage, AnalysisInfoAdmin)
admin.site.register(Residuals)
admin.site.register(ZScore)
