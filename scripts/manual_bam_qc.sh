#!/usr/bin/env sh
# manual_bam_qc.sh
# manual ATACseq bam QC

if [ $# -ne 5 ]
then
    echo "manual_bam_qc.sh sid unfiltered_bam filtered_bam peak outdir"
    exit 1
fi

sid=$1
rawbam=$2
filtbam=$3
peak=$4
outdir=$5

fragfile=${outdir}/fragment_length_count.${sid}.txt
mitofile=${outdir}/mito_content.${sid}.txt
tmpbam=${outdir}/${sid}.primary.name_sorted.bam
libfile=${outdir}/library_complexity.${sid}.txt
fripfile=${outdir}/frip.${sid}.txt

# Fragment/Insert size
samtools view ${rawbam} |awk '$9>0' |cut -f 9 |sort |uniq -c |sort -b -k2,2n |sed -e 's/^[ \t]*//' >${fragfile}

# %mitochondrial reads
mtReads=$(samtools idxstats ${rawbam} |grep 'MT' | cut -f 3)
totalReads=$(samtools idxstats ${rawbam} |awk '{SUM += $3} END {print SUM}')
mito=$(awk "BEGIN {print "${mtReads}"/"${totalReads}"}")
echo -e "mtReads\ttotalReads\tmtPercent\n${mtReads}\t${totalReads}\t${mito}" >${mitofile}

# Library complexity
# !!! the bam used here is sorted bam after duplicates marking !!! 
# remove non-primary reads, sort bam by names
samtools view -@ 8 -F 256 -h -b -O BAM ${rawbam} |samtools sort -@ 8 -n -O BAM -o ${tmpbam} -T ${outdir}/tmp -

echo -e "Total\tDistinct\tM1\tM2\tNRF\tPBC1\tPBC2" >${libfile}
bedtools bamtobed -bedpe -i ${tmpbam} | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'MT' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0}($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n", mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' >>${libfile}

# Fraction of reads in peaks (FRiP)
total_reads=$(samtools view -c ${filtbam})

reads_in_peak=$(bedtools sort -i ${peak} |bedtools merge -i stdin |bedtools intersect -u -nonamecheck -a ${filtbam} -b stdin -ubam |samtools view -c)

FRiP=$(awk "BEGIN {print "${reads_in_peak}"/"${total_reads}"}")

echo -e "readsInPeak\ttotalReads\tFRiP\n${reads_in_peak}\t${total_reads}\t${FRiP}" >${fripfile}

