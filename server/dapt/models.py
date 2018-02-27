from django.db import models


# Create your models here.

class Run(models.Model):
    pool_number = models.IntegerField(verbose_name='仪器Pooling编号')
    pool_lab = models.IntegerField(verbose_name='实验室Pooling编号')
    home_dir = models.CharField(max_length=200, verbose_name='目录全名')

    def __str__(self):
        return self.home_dir


class Sample(models.Model):
    run = models.ForeignKey(Run, verbose_name='Pool', null=True, on_delete=models.SET_NULL)
    name = models.CharField(max_length=100, verbose_name='样本编号')
    seq_barcode = models.CharField(max_length=100, verbose_name='上机条码号')
    bam_file = models.CharField(max_length=500, verbose_name='原始bam文件')
    barcode = models.CharField(max_length=100, blank=True, null=True, verbose_name='样本条码号')
    is_control = models.BooleanField(default=False, verbose_name='是否为对照样本')
    is_abnormal = models.BooleanField(default=False, verbose_name='异常样本')

    def __str__(self):
        return self.name


class ClinicalInfo(models.Model):
    sample = models.OneToOneField(Sample, verbose_name='样本', on_delete=models.CASCADE)
    age = models.IntegerField(verbose_name='年龄')
    pregnant_week = models.CharField(max_length=50, verbose_name='孕周')
    BMI = models.IntegerField(verbose_name='BMI')
    baby_num = models.TextField(max_length=50, verbose_name='胎数')
    NT = models.TextField(max_length=20, verbose_name='NT值')
    T21 = models.TextField(max_length=20, verbose_name='T21初筛')
    T18 = models.TextField(max_length=20, verbose_name='T18初筛')
    abnormal_pregnant = models.TextField(verbose_name='不良孕史')
    tube_baby = models.CharField(max_length=20, verbose_name='试管婴儿')
    all_tran = models.CharField(max_length=20, verbose_name='异体输血')
    transplantation = models.TextField(verbose_name='移植手术')
    tumor = models.TextField(verbose_name='肿瘤')


class AnalysisInfo(models.Model):
    sample = models.OneToOneField(Sample, verbose_name='样本', on_delete=models.CASCADE)
    gc_content = models.FloatField(verbose_name='gc含量(%)')
    usable_reads = models.IntegerField(verbose_name='可用reads数')
    aneuploid = models.BooleanField(default=False, verbose_name='是否为非整倍体')
    T21 = models.BooleanField(default=False, verbose_name='21号染色体异常')
    T18 = models.BooleanField(default=False, verbose_name='18号染色体异常')
    T13 = models.BooleanField(default=False, verbose_name='13号染色体异常')

    def __str__(self):
        return self.sample.name


class GCContent(models.Model):
    sample = models.OneToOneField(Sample, verbose_name='样本', on_delete=models.CASCADE)
    chr1 = models.FloatField(verbose_name='Chr1 gc content(%)')
    chr2 = models.FloatField(verbose_name='Chr2 gc content(%)')
    chr3 = models.FloatField(verbose_name='Chr3 gc content(%)')
    chr4 = models.FloatField(verbose_name='Chr4 gc content(%)')
    chr5 = models.FloatField(verbose_name='Chr5 gc content(%)')
    chr6 = models.FloatField(verbose_name='Chr6 gc content(%)')
    chr7 = models.FloatField(verbose_name='Chr7 gc content(%)')
    chr8 = models.FloatField(verbose_name='Chr8 gc content(%)')
    chr9 = models.FloatField(verbose_name='Chr9 gc content(%)')
    chr10 = models.FloatField(verbose_name='Chr10 gc content(%)')
    chr11 = models.FloatField(verbose_name='Chr11 gc content(%)')
    chr12 = models.FloatField(verbose_name='Chr12 gc content(%)')
    chr13 = models.FloatField(verbose_name='Chr13 gc content(%)')
    chr14 = models.FloatField(verbose_name='Chr14 gc content(%)')
    chr15 = models.FloatField(verbose_name='Chr15 gc content(%)')
    chr16 = models.FloatField(verbose_name='Chr16 gc content(%)')
    chr17 = models.FloatField(verbose_name='Chr17 gc content(%)')
    chr18 = models.FloatField(verbose_name='Chr18 gc content(%)')
    chr19 = models.FloatField(verbose_name='Chr19 gc content(%)')
    chr20 = models.FloatField(verbose_name='Chr20 gc content(%)')
    chr21 = models.FloatField(verbose_name='Chr21 gc content(%)')
    chr22 = models.FloatField(verbose_name='Chr22 gc content(%)')
    chrX = models.FloatField(verbose_name='ChrX gc content(%)')
    chrY = models.FloatField(verbose_name='ChrY gc content(%)')

    def __str__(self):
        return f'gc content of {self.sample.name}'


class RawCoverage(models.Model):
    sample = models.OneToOneField(Sample, verbose_name='样本', on_delete=models.CASCADE)
    chr1 = models.FloatField(verbose_name='Chr1 Raw Coverage(%)')
    chr2 = models.FloatField(verbose_name='Chr2 Raw Coverage(%)')
    chr3 = models.FloatField(verbose_name='Chr3 Raw Coverage(%)')
    chr4 = models.FloatField(verbose_name='Chr4 Raw Coverage(%)')
    chr5 = models.FloatField(verbose_name='Chr5 Raw Coverage(%)')
    chr6 = models.FloatField(verbose_name='Chr6 Raw Coverage(%)')
    chr7 = models.FloatField(verbose_name='Chr7 Raw Coverage(%)')
    chr8 = models.FloatField(verbose_name='Chr8 Raw Coverage(%)')
    chr9 = models.FloatField(verbose_name='Chr9 Raw Coverage(%)')
    chr10 = models.FloatField(verbose_name='Chr10 Raw Coverage(%)')
    chr11 = models.FloatField(verbose_name='Chr11 Raw Coverage(%)')
    chr12 = models.FloatField(verbose_name='Chr12 Raw Coverage(%)')
    chr13 = models.FloatField(verbose_name='Chr13 Raw Coverage(%)')
    chr14 = models.FloatField(verbose_name='Chr14 Raw Coverage(%)')
    chr15 = models.FloatField(verbose_name='Chr15 Raw Coverage(%)')
    chr16 = models.FloatField(verbose_name='Chr16 Raw Coverage(%)')
    chr17 = models.FloatField(verbose_name='Chr17 Raw Coverage(%)')
    chr18 = models.FloatField(verbose_name='Chr18 Raw Coverage(%)')
    chr19 = models.FloatField(verbose_name='Chr19 Raw Coverage(%)')
    chr20 = models.FloatField(verbose_name='Chr20 Raw Coverage(%)')
    chr21 = models.FloatField(verbose_name='Chr21 Raw Coverage(%)')
    chr22 = models.FloatField(verbose_name='Chr22 Raw Coverage(%)')
    chrX = models.FloatField(verbose_name='ChrX Raw Coverage(%)')
    chrY = models.FloatField(verbose_name='ChrY Raw Coverage(%)')

    def __str__(self):
        return f'Raw Coverage of {self.sample.name}'


class RawCoverage2(models.Model):
    sample = models.OneToOneField(Sample, verbose_name='样本', on_delete=models.CASCADE)
    chr1 = models.FloatField(verbose_name='Chr1 Raw Coverage(%)')
    chr2 = models.FloatField(verbose_name='Chr2 Raw Coverage(%)')
    chr3 = models.FloatField(verbose_name='Chr3 Raw Coverage(%)')
    chr4 = models.FloatField(verbose_name='Chr4 Raw Coverage(%)')
    chr5 = models.FloatField(verbose_name='Chr5 Raw Coverage(%)')
    chr6 = models.FloatField(verbose_name='Chr6 Raw Coverage(%)')
    chr7 = models.FloatField(verbose_name='Chr7 Raw Coverage(%)')
    chr8 = models.FloatField(verbose_name='Chr8 Raw Coverage(%)')
    chr9 = models.FloatField(verbose_name='Chr9 Raw Coverage(%)')
    chr10 = models.FloatField(verbose_name='Chr10 Raw Coverage(%)')
    chr11 = models.FloatField(verbose_name='Chr11 Raw Coverage(%)')
    chr12 = models.FloatField(verbose_name='Chr12 Raw Coverage(%)')
    chr13 = models.FloatField(verbose_name='Chr13 Raw Coverage(%)')
    chr14 = models.FloatField(verbose_name='Chr14 Raw Coverage(%)')
    chr15 = models.FloatField(verbose_name='Chr15 Raw Coverage(%)')
    chr16 = models.FloatField(verbose_name='Chr16 Raw Coverage(%)')
    chr17 = models.FloatField(verbose_name='Chr17 Raw Coverage(%)')
    chr18 = models.FloatField(verbose_name='Chr18 Raw Coverage(%)')
    chr19 = models.FloatField(verbose_name='Chr19 Raw Coverage(%)')
    chr20 = models.FloatField(verbose_name='Chr20 Raw Coverage(%)')
    chr21 = models.FloatField(verbose_name='Chr21 Raw Coverage(%)')
    chr22 = models.FloatField(verbose_name='Chr22 Raw Coverage(%)')
    chrX = models.FloatField(verbose_name='ChrX Raw Coverage(%)')
    chrY = models.FloatField(verbose_name='ChrY Raw Coverage(%)')

    def __str__(self):
        return f'Raw Coverage of {self.sample.name}'


class FitCoverage(models.Model):
    sample = models.OneToOneField(Sample, verbose_name='样本', on_delete=models.CASCADE)
    chr1 = models.FloatField(verbose_name='Chr1 Coverage(%)')
    chr2 = models.FloatField(verbose_name='Chr2 Coverage(%)')
    chr3 = models.FloatField(verbose_name='Chr3 Coverage(%)')
    chr4 = models.FloatField(verbose_name='Chr4 Coverage(%)')
    chr5 = models.FloatField(verbose_name='Chr5 Coverage(%)')
    chr6 = models.FloatField(verbose_name='Chr6 Coverage(%)')
    chr7 = models.FloatField(verbose_name='Chr7 Coverage(%)')
    chr8 = models.FloatField(verbose_name='Chr8 Coverage(%)')
    chr9 = models.FloatField(verbose_name='Chr9 Coverage(%)')
    chr10 = models.FloatField(verbose_name='Chr10 Coverage(%)')
    chr11 = models.FloatField(verbose_name='Chr11 Coverage(%)')
    chr12 = models.FloatField(verbose_name='Chr12 Coverage(%)')
    chr13 = models.FloatField(verbose_name='Chr13 Coverage(%)')
    chr14 = models.FloatField(verbose_name='Chr14 Coverage(%)')
    chr15 = models.FloatField(verbose_name='Chr15 Coverage(%)')
    chr16 = models.FloatField(verbose_name='Chr16 Coverage(%)')
    chr17 = models.FloatField(verbose_name='Chr17 Coverage(%)')
    chr18 = models.FloatField(verbose_name='Chr18 Coverage(%)')
    chr19 = models.FloatField(verbose_name='Chr19 Coverage(%)')
    chr20 = models.FloatField(verbose_name='Chr20 Coverage(%)')
    chr21 = models.FloatField(verbose_name='Chr21 Coverage(%)')
    chr22 = models.FloatField(verbose_name='Chr22 Coverage(%)')
    chrX = models.FloatField(verbose_name='ChrX Coverage(%)')
    chrY = models.FloatField(verbose_name='ChrY Coverage(%)')

    def __str__(self):
        return f'Coverage of {self.sample.name}'


class FitCoverage2(models.Model):
    """
    fit coverage with gc contents for every chromosome
    """
    sample = models.OneToOneField(Sample, verbose_name='样本', on_delete=models.CASCADE)
    chr1 = models.FloatField(verbose_name='Chr1 Coverage(%)')
    chr2 = models.FloatField(verbose_name='Chr2 Coverage(%)')
    chr3 = models.FloatField(verbose_name='Chr3 Coverage(%)')
    chr4 = models.FloatField(verbose_name='Chr4 Coverage(%)')
    chr5 = models.FloatField(verbose_name='Chr5 Coverage(%)')
    chr6 = models.FloatField(verbose_name='Chr6 Coverage(%)')
    chr7 = models.FloatField(verbose_name='Chr7 Coverage(%)')
    chr8 = models.FloatField(verbose_name='Chr8 Coverage(%)')
    chr9 = models.FloatField(verbose_name='Chr9 Coverage(%)')
    chr10 = models.FloatField(verbose_name='Chr10 Coverage(%)')
    chr11 = models.FloatField(verbose_name='Chr11 Coverage(%)')
    chr12 = models.FloatField(verbose_name='Chr12 Coverage(%)')
    chr13 = models.FloatField(verbose_name='Chr13 Coverage(%)')
    chr14 = models.FloatField(verbose_name='Chr14 Coverage(%)')
    chr15 = models.FloatField(verbose_name='Chr15 Coverage(%)')
    chr16 = models.FloatField(verbose_name='Chr16 Coverage(%)')
    chr17 = models.FloatField(verbose_name='Chr17 Coverage(%)')
    chr18 = models.FloatField(verbose_name='Chr18 Coverage(%)')
    chr19 = models.FloatField(verbose_name='Chr19 Coverage(%)')
    chr20 = models.FloatField(verbose_name='Chr20 Coverage(%)')
    chr21 = models.FloatField(verbose_name='Chr21 Coverage(%)')
    chr22 = models.FloatField(verbose_name='Chr22 Coverage(%)')
    chrX = models.FloatField(verbose_name='ChrX Coverage(%)')
    chrY = models.FloatField(verbose_name='ChrY Coverage(%)')

    def __str__(self):
        return f'Coverage of {self.sample.name}'


class Residuals(models.Model):
    sample = models.OneToOneField(Sample, verbose_name='样本', on_delete=models.CASCADE)
    chr1 = models.FloatField(verbose_name='Chr1 Residuals')
    chr2 = models.FloatField(verbose_name='Chr2 Residuals')
    chr3 = models.FloatField(verbose_name='Chr3 Residuals')
    chr4 = models.FloatField(verbose_name='Chr4 Residuals')
    chr5 = models.FloatField(verbose_name='Chr5 Residuals')
    chr6 = models.FloatField(verbose_name='Chr6 Residuals')
    chr7 = models.FloatField(verbose_name='Chr7 Residuals')
    chr8 = models.FloatField(verbose_name='Chr8 Residuals')
    chr9 = models.FloatField(verbose_name='Chr9 Residuals')
    chr10 = models.FloatField(verbose_name='Chr10 Residuals')
    chr11 = models.FloatField(verbose_name='Chr11 Residuals')
    chr12 = models.FloatField(verbose_name='Chr12 Residuals')
    chr13 = models.FloatField(verbose_name='Chr13 Residuals')
    chr14 = models.FloatField(verbose_name='Chr14 Residuals')
    chr15 = models.FloatField(verbose_name='Chr15 Residuals')
    chr16 = models.FloatField(verbose_name='Chr16 Residuals')
    chr17 = models.FloatField(verbose_name='Chr17 Residuals')
    chr18 = models.FloatField(verbose_name='Chr18 Residuals')
    chr19 = models.FloatField(verbose_name='Chr19 Residuals')
    chr20 = models.FloatField(verbose_name='Chr20 Residuals')
    chr21 = models.FloatField(verbose_name='Chr21 Residuals')
    chr22 = models.FloatField(verbose_name='Chr22 Residuals')
    chrX = models.FloatField(verbose_name='ChrX Residuals')
    chrY = models.FloatField(verbose_name='ChrY Residuals')

    def __str__(self):
        return f'Residuals of {self.sample.name}'


class ZScore(models.Model):
    sample = models.OneToOneField(Sample, verbose_name='样本', on_delete=models.CASCADE)
    chr1 = models.FloatField(verbose_name='Chr1 Z-Score')
    chr2 = models.FloatField(verbose_name='Chr2 Z-Score')
    chr3 = models.FloatField(verbose_name='Chr3 Z-Score')
    chr4 = models.FloatField(verbose_name='Chr4 Z-Score')
    chr5 = models.FloatField(verbose_name='Chr5 Z-Score')
    chr6 = models.FloatField(verbose_name='Chr6 Z-Score')
    chr7 = models.FloatField(verbose_name='Chr7 Z-Score')
    chr8 = models.FloatField(verbose_name='Chr8 Z-Score')
    chr9 = models.FloatField(verbose_name='Chr9 Z-Score')
    chr10 = models.FloatField(verbose_name='Chr10 Z-Score')
    chr11 = models.FloatField(verbose_name='Chr11 Z-Score')
    chr12 = models.FloatField(verbose_name='Chr12 Z-Score')
    chr13 = models.FloatField(verbose_name='Chr13 Z-Score')
    chr14 = models.FloatField(verbose_name='Chr14 Z-Score')
    chr15 = models.FloatField(verbose_name='Chr15 Z-Score')
    chr16 = models.FloatField(verbose_name='Chr16 Z-Score')
    chr17 = models.FloatField(verbose_name='Chr17 Z-Score')
    chr18 = models.FloatField(verbose_name='Chr18 Z-Score')
    chr19 = models.FloatField(verbose_name='Chr19 Z-Score')
    chr20 = models.FloatField(verbose_name='Chr20 Z-Score')
    chr21 = models.FloatField(verbose_name='Chr21 Z-Score')
    chr22 = models.FloatField(verbose_name='Chr22 Z-Score')
    chrX = models.FloatField(verbose_name='ChrX Z-Score')
    chrY = models.FloatField(verbose_name='ChrY Z-Score')

    def __str__(self):
        return f'Z-Score of {self.sample.name}'
