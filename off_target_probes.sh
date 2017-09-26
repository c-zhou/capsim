#!/bin/bash

printUsage() {
    echo ""
    echo ""
    echo "  Usage: bash off_target_probes.sh [options] -b <bed> -r <fasta> -a <bam> -q <fastq>" 
    echo "  Options:"
    echo "      -b/--target-bed                  Bed file of the targe regions. "
    echo "      -r/--reference                   Reference genome fasta file. "
    echo "      -a/--bam                         XXX. "
    echo "      -w/--window-size                 Window size for statistics of the depth of coverage of"
    echo "                                       the off targe regions (default 1000)."
    echo "      -d/--min-depth                   Minimum depth of coverage of the off target regions to"
    echo "                                       analyse (default 10000). "
    echo "      -q/--probe-seq                   The fastq file for probe sequences."
    echo "      -t/--threads                     Number of threads for alignment (default 1)."
    echo "      -p/--prefix                      Prefix of the output files (default ./out)."
    echo ""
}

TARGET_BED=""
REFERENCE=""
BAM=""
WINDOW_SIZE=1000
MIN_DEPTH=10000
PROBE_SEQ=""
THREADS=1
PREFIX_OUT="./out"

while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
        -b|--target-bed)
            TARGET_BED="$2"
            shift
            shift
            ;;
        -r|--reference)
            REFERENCE="$2"
            shift
            shift
            ;;
        -a|--bam)
            BAM="$2"
            shift
            shift
            ;;
        -w|--window-size)
            WINDOW_SIZE="$2"
            shift
            shift
            ;;
        -d|--min-depth)
            MIN_DEPTH="$2"
            shift
            shift
            ;;
        -q|--probe_seq)
            PROBE_SEQ="$2"
            shift
            shift
            ;;
        -t|--threads)
            THREADS="$2"
            shift
            shift
            ;;
        -p|--prefix)
            PREFIX_OUT="$2"
            shift
            shift
            ;;
        *) # unknown option
            echo ""
            echo "!!!ERROR: unknown options $1"
            printUsage
            exit
            ;;
    esac
done

if [[ $TARGET_BED == "" ]]
then
    echo ""
    echo "!!!ERROR: A bed file for targe regions is required: -b/--target-bed"
    printUsage
    exit
fi

if [[ $REFERENCE == "" ]]
then
    echo ""
    echo "!!!ERROR: A reference genome fasta file is required: -r/--reference"
    printUsage
    exit
fi

if [[ $BAM == "" ]]
then
    echo ""
    echo "!!!ERROR: A bam file XXX is required: -a/--bam"
    printUsage
    exit
fi

if [[ $PROBE_SEQ == "" ]]
then
    echo ""
    echo "!!!ERROR: A fastq file for probe sequences is required: -q/--probe-seq"
    printUsage
    exit
fi

echo ""
echo "  Options provided:"
echo "      TARGET BED     = $TARGET_BED"
echo "      REFERENCE      = $REFERENCE"
echo "      BAM            = $BAM"
echo "      WINDOW SIZE    = $WINDOW_SIZE"
echo "      MIN DEPTH      = $MIN_DEPTH"
echo "      PROBE SEQ      = $PROBE_SEQ"
echo "      THREADS        = $THREADS"
echo "      PREFIX OUT     = $PREFIX_OUT"

: <<"SKIP"
bedtools complement -i Target_regions_500.bed -g hg19.genome > Target_regions_500_complement.bed

bedtools intersect -abam MS_hg19_primary_sorted.bam -b Target_regions_500_complement.bed > MS_hg19_primary_sorted_target_500_complement.bam

samtools index MS_hg19_primary_sorted_target_500_complement.bam MS_hg19_primary_sorted_target_500_complement.bai

awk 'BEGIN{step=1000}{num=($3-$2)/step; for(i=0;i<num-1;i++) print $1,$2+step*i,$2+step*(i+1)-1; if(step*int(num)!=$3-$2) print $1,$2+step*int(num),$3;}' Target_regions_500_complement.bed > Target_regions_500_complement_1Kb.bed 

samtools bedcov Target_regions_500_complement_1Kb.bed MS_hg19_primary_sorted_target.bam > MS_hg19_primary_sorted_target_complement_1kb.bedcov

awk '{if ($4>10000) {print $1"\t"$2"\t"$3"\t"}}' MS_hg19_primary_sorted_target_complement_1kb.bedcov > MS_off_target_regions.bed

bedtools getfasta -fi hg19.fasta -bed MS_off_target_regions.bed -fo MS_off_target_regions.fas

bwa index -a bwtsw -p MS_off_target_regions MS_off_target_regions.fas

bwa mem -t 8 -Y -c 1000 MS_off_target_regions probes.fastq | samtools view -@8 -bS | samtools sort -@8 -o MS_off-target_alignement_sorted.bam && samtools index MS_off-target_alignement_sorted.bam MS_off-target_alignement_sorted.bai

samtools view -F4 MS_off-target_alignement_sorted.bam | cut -f1 | sort -u > MS_off-target_alignement_probes.txt

SKIP
