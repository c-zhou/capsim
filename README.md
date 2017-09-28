# capsim

Capsim is a tool to simulate capture sequencing. It was developed as part
of the Japsa package. 

## Quick installation guide:
On a Linux/Mac machine with `make' and `git' installed, the software can be installed with

    git clone https://github.com/mdcao/japsa
    cd japsa
    make install \
      [INSTALL_DIR=~/.usr/local \] 
      [MXMEM=7000m \] 
      [SERVER=true \]


Details of installation (including for Windows) and usage of Japsa can be found 
in its documentation hosted on [ReadTheDocs](http://japsa.readthedocs.org/en/latest/index.html) 

## Usage

Before running capsim, the probe sequences need to be aligned to the reference sequence (or the genome sequence to simulate sequencing). We recommend using bowtie2 for the alignment.
   
    #Skip this step if the index has been generated
    bowtie2-build ref.fasta ref
    
    #Align the probes into the reference
    bowtie2 --local --very-sensitive-local --mp 8 --rdg 10,8 --rfg 10,8 -k 10000 -f -x ref -U probes.fa -S probes.sam

Note: for some reason, bowtie2 only accepts the query fasta file (probes.fa) containing one sequence per line.

After alignment, sort the and index the bam file with samtools

     samtools view -bSU probes.sam | samtools sort -o probes.bam -

Capsim takes the bam file as the input:

    jsa.sim.capsim --reference ref.fasta --probe probes.bam --ID someid --fmedian 5000  --pacbio output --pblen 3000 --num 20000000

or 

    jsa.sim.capsim --reference ref.fasta --probe probes.bam --ID someid --fmedian 500  --miseq output --illen 300 --num 20000000


## Subsequent analysis

A script file _off_target_probes.sh_ used to conduct the subsequent analysis to XXX is provided in this repository. To run this script file,

    bash off_target_probes.sh -b target_regions.bed -r ref.fasta -a ms_sorted.bam -w 1000 -d 10000 -q probes.fq -t 4 -p out

where,

<pre>
-b/--target-bed                  Bed file of the targe regions. 
-r/--reference                   Reference genome fasta file. 
-a/--bam                         XXX. 
-w/--window-size                 Window size for statistics of the depth of coverage of
                                 the off targe regions (default 1000).
-d/--min-depth                   Minimum depth of coverage of the off target regions to
                                 analyse (default 10000). 
-q/--probe-seq                   The file for probe sequences.
-t/--threads                     Number of threads for alignment (default 1).
-p/--prefix                      Prefix of the output files (default ./out).
</pre>
    
The following tools should be installed and added to system path.

* [samtools](http://samtools.sourceforge.net/)
* [bwa](http://bio-bwa.sourceforge.net/)
* [bedtools](https://github.com/arq5x/bedtools2)

**_The script file consists of 11 main commands which could be run step by step._**

1\. generate reference index file,
    
    samtools faidx ref.fasta

2\. generate whole genome bed file,

    awk '{print $1"\t"$2}' ref.fasta > out_ref.fasta.bed

3\. generate a bed file represents the complementary regions of the target regions,

    bedtools complement -i target_regions.bed -g ref.fasta > out_target_regions_complement.bed

2\. generate a bam file contains alignment records that mapped to the complementary regions,

    bedtools intersect -abam ms_sorted.bam -b out_target_regions_complement.bed > out_ms_sorted_complement.bam
    
3\. generate index file for the bam file,

    samtools index out_ms_sorted_complement.bam out_ms_sorted_complement.bai
    
4\. generate a bed file represents 1Kb windows of the complementary regions,

    awk -v w=1000 '{num=($3-$2)/w; for(i=0;i<num-1;i++) print $1"\t"($2+w*i)"\t"($2+w*(i+1)-1); if(w*int(num)!=$3-$2) print $1"\t"($2+w*int(num))"\t"$3;}' out_target_regions_complement.bed > out_target_regions_complement_1Kb.bed

This file will be used to calculate the depth of coverage of the bam file across the complemetary regions. Window size other than 1Kb could be used here. A smaller window size will generally result in more precise statistics but will be more time-consuming. The window size could be changed by the awk parameter _w_.

5\. calculate the depth of coverage of the bam file across the complemetary regions,

    samtools bedcov out_target_regions_complement_1Kb.bed out_ms_sorted_complement.bam > out_ms_sorted_complement_1Kb.bedcov
    
6\. filter the complementary regions according to the depth of coverage,

    awk -v d=10000 '{if ($4>d) print $1"\t"$2"\t"$3}' out_ms_sorted_complement_1Kb.bedcov > out_ms_off_target_regions.bed
    
The threshold of the depth of coverage could be specified by the awk parameter _d_.

7\. generate the fasta file of the off target regions,

    bedtools getfasta -fi ref.fasta -bed out_ms_off_target_regions.bed -fo out_ms_off_target_regions.fasta
    
8\. generate bwa index files for the off target region fasta file,

    bwa index -a bwtsw -p out_ms_off_target_regions out_ms_off_target_regions.fasta

9\. align the probe sequences to the off target region fasta file,

    bwa mem -t 8 -Y -c 1000 out_ms_off_target_regions probes.fastq | samtools view -@8 -bS | samtools sort -@8 -o out_ms_off_target_regions.bam && samtools index out_ms_off_target_regions.bam out_ms_off_target_regions.bai
    
10\. extract the off target alignment probe sequences,

    samtools view -F4 out_ms_off_target_regions.bam | cut -f1 | sort -u > out_ms_off_target_alignment_probes.txt
    

