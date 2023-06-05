# benchmarking-structural-variant-callers
### Description
In this project, I benchmarked several long-read structural variant callers (cuteSV, DeBreak, PBSV, Sniffles, SVDSS, and SVIM) on a range of sequencing depths (5X, 10X, 15X, 20X, 25X, 30X, MAX) in order to determine the minimum depth allowing for the most sufficient yield of structural variant calls. With these benchmarking statistics, a cost-effective sequencing depth can be determined for use in the generation of new data. For this analysis, I used data of the SK-BR-3 breast cancer cell line generated on three sequencing platforms: Oxford Nanopore Technologies (ONT), PacBio Continuous Long Reads (PBCLR), and Illumina short-reads (SR).
 ## 1. Alignment
 I used minimap2 to align long-reads (ONT and PBCLR) and bwa-mem2 to align the short-reads (SR) onto the GRCh38.p13 assembly.
 
```
# indexing

minimap2 -d GRCh38.p13.genome.mmi GRCh38.p13.genome.fa
bwa-mem2 index -p bwa_index GRCh38.p13.genome.fa

# alignment

minimap2 -ax map-ont -t 24 -Y -R '@RG\tID:SKBR3\tPL:ONT\tLB:library\tSM:ONT' \
GRCh38.p13.genome.mmi ONT.fastq > SKBR3_ONT_MAX.sam

minimap2 -ax map-pb -t 24 -Y -R '@RG\tID:SKBR3\tPL:PACBIO\tLB:library\tSM:PBCLR' \
GRCh38.p13.genome.mmi PBCLR.fastq > SKBR3_PBCLR_MAX.sam

bwa-mem2 mem -t 24 -M -Y -R '@RG\tID:SKBR3\tPL:ILLUMINA\tLB:library\tSM:SR' \
bwa_index SR_R1.fastq SR_R2.fastq > SKBR3_SR_MAX.sam
```
The -Y and -R options are both necessary for downstream SV analysis. -Y determines softclipping for supplementary alignments. -R is required by PBSV and GATK tools. [See here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups) learn about read groups. 

After the alignment, the convert the three SAM files into sorted BAMs. 
```
for samfile in *MAX.sam
do
name=$(basename ${samfile} _MAX.sam)
samtools sort ${samfile} -o ${name}_sort_MAX.bam -@ 16
done
```
And the resulting filenames would be:
>SKBR3_ONT_sort_MAX.bam<br>
>SKBR3_PBCLR_sort_MAX.bam<br>
>SKBR3_SR_sort_MAX.bam<br>


