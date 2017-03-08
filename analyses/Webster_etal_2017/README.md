# Analyses from Webster et al. 2017
This directory contains scripts and information for reproducing the analyses that
accompany the initial publication of XYalign (Webster et al., 2017).

## Data
For these analyses, we used three different publicly-available datasets:

1) Exome, low-coverage whole-genome, and high-coverage whole-genome data from
two individuals in the 1000 genomes project (HG00512 - male, HG00513 - female)

Fastqs downloaded from:
```
HG00512

Exome
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR034/­ERR034508/­ERR034508_1.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR034/­ERR034508/­ERR034508_2.­fastq.­gz

Low-coverage whole-genome
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR016/­ERR016114/­ERR016114_1.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR016/­ERR016114/­ERR016114_2.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR016/­ERR016116/­ERR016116_1.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR016/­ERR016116/­ERR016116_2.­fastq.­gz

High-coverage whole genome
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR894/­ERR894725/­ERR894725_1.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR894/­ERR894725/­ERR894725_2.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR894/­ERR894726/­ERR894726_1.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR894/­ERR894726/­ERR894726_2.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR899/­ERR899712/­ERR899712_1.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR899/­ERR899712/­ERR899712_2.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR899/­ERR899713/­ERR899713_1.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR899/­ERR899713/­ERR899713_2.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR903/­ERR903029/­ERR903029_1.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR903/­ERR903029/­ERR903029_2.­fastq.­gz

HG00513

Exome
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR034/­ERR034509/­ERR034509_1.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR034/­ERR034509/­ERR034509_2.­fastq.­gz

Low-coverage whole-genome
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR016/­ERR016118/­ERR016118_1.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR016/­ERR016118/­ERR016118_2.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR016/­ERR016119/­ERR016119_1.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR016/­ERR016119/­ERR016119_2.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR016/­ERR016121/­ERR016121_1.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR016/­ERR016121/­ERR016121_2.­fastq.­gz

High-coverage whole-genome
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR894/­ERR894727/­ERR894727_1.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR894/­ERR894727/­ERR894727_2.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR894/­ERR894728/­ERR894728_1.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR894/­ERR894728/­ERR894728_2.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR899/­ERR899714/­ERR899714_1.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR899/­ERR899714/­ERR899714_2.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR899/­ERR899715/­ERR899715_2.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR899/­ERR899715/­ERR899715_1.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR899/­ERR899716/­ERR899716_1.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR899/­ERR899716/­ERR899716_2.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR903/­ERR903027/­ERR903027_1.­fastq.­gz
ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR903/­ERR903027/­ERR903027_2.­fastq.­gz
```

2) 24 high-coverage whole-genomes from the 1000 genomes project

BAM files downloaded from:

```
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/high_coverage_alignment/HG00096.wgs.ILLUMINA.bwa.GBR.high_cov_pcr_free.20140203.bam
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG01879/high_coverage_alignment/HG01879.wgs.ILLUMINA.bwa.ACB.high_cov_pcr_free.20140203.bam
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG03052/high_coverage_alignment/HG03052.wgs.ILLUMINA.bwa.MSL.high_cov_pcr_free.20140203.bam
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG01051/high_coverage_alignment/HG01051.wgs.ILLUMINA.bwa.PUR.high_cov_pcr_free.20140203.bam
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19625/high_coverage_alignment/NA19625.wgs.ILLUMINA.bwa.ASW.high_cov_pcr_free.20140203.bam
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG01583/high_coverage_alignment/HG01583.wgs.ILLUMINA.bwa.PJL.high_cov_pcr_free.20140203.bam
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG03742/high_coverage_alignment/HG03742.wgs.ILLUMINA.bwa.ITU.high_cov_pcr_free.20140203.bam
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA20502/high_coverage_alignment/NA20502.wgs.ILLUMINA.bwa.TSI.high_cov_pcr_free.20140203.bam
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18525/high_coverage_alignment/NA18525.wgs.ILLUMINA.bwa.CHB.high_cov_pcr_free.20140203.bam
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG02922/high_coverage_alignment/HG02922.wgs.ILLUMINA.bwa.ESN.high_cov_pcr_free.20140203.bam
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19648/high_coverage_alignment/NA19648.wgs.ILLUMINA.bwa.MXL.high_cov_pcr_free.20140203.bam
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00419/high_coverage_alignment/HG00419.wgs.ILLUMINA.bwa.CHS.high_cov_pcr_free.20140203.bam
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG01112/high_coverage_alignment/HG01112.wgs.ILLUMINA.bwa.CLM.high_cov_pcr_free.20140203.bam
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19017/high_coverage_alignment/NA19017.wgs.ILLUMINA.bwa.LWK.high_cov_pcr_free.20140203.bam
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00268/high_coverage_alignment/HG00268.wgs.ILLUMINA.bwa.FIN.high_cov_pcr_free.20140203.bam
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG02568/high_coverage_alignment/HG02568.wgs.ILLUMINA.bwa.GWD.high_cov_pcr_free.20140203.bam
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG01500/high_coverage_alignment/HG01500.wgs.ILLUMINA.bwa.IBS.high_cov_pcr_free.20140203.bam
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18939/high_coverage_alignment/NA18939.wgs.ILLUMINA.bwa.JPT.high_cov_pcr_free.20140203.bam
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG03642/high_coverage_alignment/HG03642.wgs.ILLUMINA.bwa.STU.high_cov_pcr_free.20140203.bam
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG03006/high_coverage_alignment/HG03006.wgs.ILLUMINA.bwa.BEB.high_cov_pcr_free.20140203.bam
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00759/high_coverage_alignment/HG00759.wgs.ILLUMINA.bwa.CDX.high_cov_pcr_free.20140203.bam
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA20845/high_coverage_alignment/NA20845.wgs.ILLUMINA.bwa.GIH.high_cov_pcr_free.20140203.bam
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG01595/high_coverage_alignment/HG01595.wgs.ILLUMINA.bwa.KHV.high_cov_pcr_free.20140203.bam
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG01565/high_coverage_alignment/HG01565.wgs.ILLUMINA.bwa.PEL.high_cov_pcr_free.20140203.bam
```

3) Six gorillas (3 males and 3 females) from two different species

Fastqs downloaded from SRA with the following command lines using [fastq-dump](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump):

Azizi - Male - Gorilla gorilla gorilla
```
for i in SRR747941 SRR747942 SRR747939 SRR747940; do fastq-dump --gzip --readids --split-files $i; done
```

Banjo - Male - Gorilla gorilla gorilla
```
for i in SRR747969 SRR747970 SRR747971 SRR747967 SRR747968; do fastq-dump --gzip --readids --split-files $i; done
```

Delphi - Female - Gorilla gorilla gorilla
```
for i in SRR747987 SRR747988 SRR747989 SRR747984 SRR747985 SRR747986; do fastq-dump --gzip --readids --split-files $i; done
```

Nyango - Female - Gorilla gorilla diehli
```
for i in SRR748109 SRR748110 SRR748111 SRR748112; do fastq-dump --gzip --readids --split-files $i; done
```
Kaisi - Male - Gorilla beringei graueri
```
for i in SRR747658 SRR747657 SRR747654 SRR747655 SRR747656 SRR747651 SRR747652 SRR747653; do fastq-dump --gzip --readids --split-files $i; done
```

Victoria - Female - Gorilla beringei graueri
```
for i in SRR748192 SRR748191 SRR748190 SRR748189 SRR748188 SRR748187; do fastq-dump --gzip --readids --split-files $i; done
```

## Directory Structure


## Running snakemake


## Citations

####XYalign
Webster TH et al. In prep. XYalign: inferring sex chromosome content and correcting
for technical biases in next-generation sequencing data.

####Data sources
