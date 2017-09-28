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
    echo "      -q/--probe-seq                   The file for probe sequences."
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

STEP=1

echo "["`date`"]" STEP $((STEP++)): Generate reference index
samtools faidx $REFERENCE 

echo "["`date`"]" STEP $((STEP++)):
awk '{print $1"\t"$2}' "$REFERENCE".fai > "$PREFIX_OUT"_$(basename $REFERENCE).genome

echo "["`date`"]" STEP $((STEP++)):
bedtools complement -i $TARGET_BED -g "$PREFIX_OUT"_$(basename $REFERENCE).genome | awk -v w="$WINDOW_SIZE" '{num=($3-$2)/w; for(i=0;i<num-1;i++) print $1"\t"($2+w*i)"\t"($2+w*(i+1)-1); if(w*int(num)!=$3-$2) print $1"\t"($2+w*int(num))"\t"$3;}' > "$PREFIX_OUT"_target_regions_complement.bed

echo "["`date`"]" STEP $((STEP++)):
bedtools intersect -abam $BAM -b "$PREFIX_OUT"_target_regions_complement.bed > "$PREFIX_OUT"_`echo $(basename $BAM) | sed 's/.bam$//g'`_target_complement.bam

echo "["`date`"]" STEP $((STEP++)):
samtools index "$PREFIX_OUT"_`echo $(basename $BAM) | sed 's/.bam$//g'`_target_complement.bam "$PREFIX_OUT"_`echo $(basename $BAM) | sed 's/.bam$//g'`_target_complement.bai

echo "["`date`"]" STEP $((STEP++)):
samtools bedcov "$PREFIX_OUT"_target_regions_complement.bed "$PREFIX_OUT"_`echo $(basename $BAM) | sed 's/.bam$//g'`_target_complement.bam | awk -v d="$MIN_DEPTH" '{if ($4>d) {print $1"\t"$2"\t"$3"\t"}}' > "$PREFIX_OUT"_`echo $(basename $BAM) | sed 's/.bam$//g'`_off_target_regions.bed

echo "["`date`"]" STEP $((STEP++)):
bedtools getfasta -fi $REFERENCE -bed "$PREFIX_OUT"_`echo $(basename $BAM) | sed 's/.bam$//g'`_off_target_regions.bed -fo "$PREFIX_OUT"_`echo $(basename $BAM) | sed 's/.bam$//g'`_off_target_regions.fas

echo "["`date`"]" STEP $((STEP++)):
bwa index -a bwtsw "$PREFIX_OUT"_`echo $(basename $BAM) | sed 's/.bam$//g'`_off_target_regions.fas 

echo "["`date`"]" STEP $((STEP++)):
tail -n +2 $PROBE_SEQ | awk '{printf(">%s:%s\n%s\n",$1,$2,$3)}' > "$PREFIX_OUT"_"$(basename $PROBE_SEQ)".fa

echo "["`date`"]" STEP $((STEP++)):
bwa mem -t $THREADS -Y -c 1000 "$PREFIX_OUT"_`echo $(basename $BAM) | sed 's/.bam$//g'`_off_target_regions.fas "$PREFIX_OUT"_"$(basename $PROBE_SEQ)".fa | samtools view -@"$THREADS" -bS | samtools sort -@"$THREADS" -o "$PREFIX_OUT"_`echo $(basename $BAM) | sed 's/.bam$//g'`_off_target_regions.bam && samtools index "$PREFIX_OUT"_`echo $(basename $BAM) | sed 's/.bam$//g'`_off_target_regions.bam "$PREFIX_OUT"_`echo $(basename $BAM) | sed 's/.bam$//g'`_off_target_regions.bai

echo "["`date`"]" STEP $((STEP++)):
samtools view -F4 "$PREFIX_OUT"_`echo $(basename $BAM) | sed 's/.bam$//g'`_off_target_regions.bam | cut -f1 | sort -u > "$PREFIX_OUT"_off_target_probes.txt

