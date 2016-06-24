### Script Collection

1. Quality filtering with Prinseq
2. Transcriptome assembly with Trinity
3. Read Mapping with BWA 

### (1) Quality Filerting with Prinseq
```HPC
soon
```

### (2) De-novo Transcriptome Assembly with Trinity

```HPC
### =================================
### Project Goby
### Scripts for sciCORE
### trinity.sh
### 27/04/16
### =================================

#!/bin/bash
#$ -N trinity_denovo
#$ -S /bin/bash
#$ -pe smp 16
#$ -l membycore=14G
#$ -l runtime=120:00:00
#$ -o /scicore/home/holmp/GROUP/01_Goby_RNAseq/z_log/trinity_qsub_across.stdout
#$ -e /scicore/home/holmp/GROUP/01_Goby_RNAseq/z_log/trinity_qsub_across.error
#$ -m beas -M jean-claude.walser@unibas.ch

## --------------------------
## load required modules
## --------------------------

module load PRINSEQ/0.20.4-goolf-1.4.10-Perl-5.16.3 Trinity/2.1.1-goolf-1.4.10 EMBOSS/6.6.0-goolf-1.4.10 

## --------------------------
## set paths
## --------------------------

gd=/scicore/home/holmp/GROUP/01_Goby_RNAseq
N50=/scicore/home/holmp/walser/bin/N50Stat.pl

mkdir -p $gd/c_assembly/trinity
date "+START Script: %T %d/%m/%Y"

## ------------------------------
## Data sets
## ------------------------------

S141R1=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/A1_141_qf_1.fastq
S141R2=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/A1_141_qf_2.fastq

S123R1=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/A2_123_qf_1.fastq
S123R2=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/A2_123_qf_2.fastq

S002R1=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/B1_2_qf_1.fastq
S002R2=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/B1_2_qf_2.fastq

S127R1=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/B2_127_qf_1.fastq
S127R2=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/B2_127_qf_2.fastq

S028R1=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/C1_28_qf_1.fastq
S028R2=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/C1_28_qf_2.fastq

S132R1=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/C2_132_qf_1.fastq
S132R2=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/C2_132_qf_2.fastq

S030R1=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/D1_30_qf_1.fastq
S030R2=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/D1_30_qf_2.fastq

S135R1=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/D2_135_qf_1.fastq
S135R2=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/D2_135_qf_2.fastq

S050R1=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/E1_50_qf_1.fastq
S050R2=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/E1_50_qf_2.fastq

S147R1=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/E2_147_qf_1.fastq
S147R2=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/E2_147_qf_2.fastq

S070R1=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/F1_70_qf_1.fastq
S070R2=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/F1_70_qf_2.fastq

S148R1=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/F2_148_qf_1.fastq
S148R2=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/F2_148_qf_2.fastq

S080R1=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/G1_80_qf_1.fastq
S080R2=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/G1_80_qf_2.fastq

S167R1=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/G2_167_qf_1.fastq
S167R2=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/G2_167_qf_2.fastq

S098R1=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/H1_98_qf_1.fastq
S098R2=/scicore/home/holmp/GROUP/01_Goby_RNAseq/b_qf/H1_98_qf_2.fastq

## -------------------------------
## Combine reads
## -------------------------------

date "+START Combine Reads: %T %d/%m/%Y"

cat $S141R1 $S123R1 $S002R1 $S127R1 $S028R1 $S132R1 $S030R1 $S135R1 $S050R1 $S147R1 $S070R1 $S148R1 $S080R1 $S167R1 $S098R1 > $TMPDIR/R1.fq
cat $S141R2 $S123R2 $S002R2 $S127R2 $S028R2 $S132R2 $S030R2 $S135R2 $S050R2 $S147R2 $S070R2 $S148R2 $S080R2 $S167R2 $S098R2 > $TMPDIR/R2.fq

## -------------------------------
## Dereplicate reads
## -------------------------------

date "+START Dereplicate Reads: %T %d/%m/%Y"

prinseq-lite.pl -verbose -out_format 3 -derep 14 -fastq $TMPDIR/R1.fq -fastq2 $TMPDIR/R2.fq -out_good $TMPDIR/R2_derep -out_bad null -log $gd/c_assembly/trinity

## -------------------------------
## De-novo Transcriptome Assembly
## -------------------------------

date "+START Trinity: %T %d/%m/%Y"

Trinity --seqType fq --min_contig_length 100 --min_glue 4 --group_pairs_distance 300 --path_reinforcement_distance 85 --min_kmer_cov 4 --max_memory 50G --CPU 16 --bflyCalculateCPU --left $TMPDIR/R2_derep.1.fastq --right $TMPDIR/R2_derep.2.fastq --output $gd/c_assembly/trinity | tee $gd/c_assembly/trinity.log

## -------------------------------
## Assembly Summary
## -------------------------------

date "+START Summary: %T %d/%m/%Y"

## Summary
$N50 -i $gd/c_assembly/trinity/Trinity.fasta -o $gd/c_assembly/trinity/Trinity.stats
## Length and %GC report for R import
infoseq -only -length -pgc $gd/c_assembly/trinity/Trinity.fasta | sed '1d' > $gd/c_assembly/trinity/Trinity_l_gc.txt

date "+End Script: %T %d/%m/%Y"
```

### (3) Read Mapping with BWA
```HPC
soon
```
