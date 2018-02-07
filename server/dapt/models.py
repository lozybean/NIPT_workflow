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
    bam_file = models.CharField(max_length=100, verbose_name='原始bam文件')
    barcode = models.CharField(max_length=100, blank=True, null=True, verbose_name='样本条码号')

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


class Coverage(models.Model):
    chrom = models.TextField(max_length=20, verbose_name='染色体')
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


class Epsilon(models.Model):
    chrom = models.TextField(max_length=20, verbose_name='染色体')
    sample = models.OneToOneField(Sample, verbose_name='样本', on_delete=models.CASCADE)
    chr1 = models.FloatField(verbose_name='Chr1 Epsilon')
    chr2 = models.FloatField(verbose_name='Chr2 Epsilon')
    chr3 = models.FloatField(verbose_name='Chr3 Epsilon')
    chr4 = models.FloatField(verbose_name='Chr4 Epsilon')
    chr5 = models.FloatField(verbose_name='Chr5 Epsilon')
    chr6 = models.FloatField(verbose_name='Chr6 Epsilon')
    chr7 = models.FloatField(verbose_name='Chr7 Epsilon')
    chr8 = models.FloatField(verbose_name='Chr8 Epsilon')
    chr9 = models.FloatField(verbose_name='Chr9 Epsilon')
    chr10 = models.FloatField(verbose_name='Chr10 Epsilon')
    chr11 = models.FloatField(verbose_name='Chr11 Epsilon')
    chr12 = models.FloatField(verbose_name='Chr12 Epsilon')
    chr13 = models.FloatField(verbose_name='Chr13 Epsilon')
    chr14 = models.FloatField(verbose_name='Chr14 Epsilon')
    chr15 = models.FloatField(verbose_name='Chr15 Epsilon')
    chr16 = models.FloatField(verbose_name='Chr16 Epsilon')
    chr17 = models.FloatField(verbose_name='Chr17 Epsilon')
    chr18 = models.FloatField(verbose_name='Chr18 Epsilon')
    chr19 = models.FloatField(verbose_name='Chr19 Epsilon')
    chr20 = models.FloatField(verbose_name='Chr20 Epsilon')
    chr21 = models.FloatField(verbose_name='Chr21 Epsilon')
    chr22 = models.FloatField(verbose_name='Chr22 Epsilon')
    chrX = models.FloatField(verbose_name='ChrX Epsilon')
    chrY = models.FloatField(verbose_name='ChrY Epsilon')

    def __str__(self):
        return f'Epsilon of {self.sample.name}'


class ZScore(models.Model):
    chrom = models.TextField(max_length=20, verbose_name='染色体')
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