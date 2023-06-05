# benchmarking-structural-variant-callers
### Description
In this project, I benchmarked several long-read structural variant callers (cuteSV, DeBreak, PBSV, Sniffles, SVDSS, and SVIM) on a range of sequencing depths (5X, 10X, 15X, 20X, 25X, 30X, MAX) in order to determine the minimum depth allowing for the most sufficient yield of structural variant calls. With these benchmarking statistics, a cost-effective sequencing depth can be determined to generate new data. For this analysis, I used data of the SK-BR-3 breast cancer cell from three sequencing platforms: Oxford Nanopore Technologies (ONT), PacBio Continuous Long Reads (PBCLR), and Illumina short-reads (SR).
 ## 1. Alignment
 I used minimap2 to align long-reads (ONT and PBCLR) and bwa-mem2 to align the short-reads (SR) using the GRCh38.p13 assembly.
```
#!/bin/bash
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

After the alignment, samtools can convert the three SAM files into sorted BAMs. 
```
for samfile in *MAX.sam
do
name=$(basename ${samfile} _MAX.sam)
samtools sort ${samfile} -o ${name}_sort_MAX.bam -@ 16
done
```
And the resulting file names would be:
>SKBR3_ONT_sort_MAX.bam<br>
>SKBR3_PBCLR_sort_MAX.bam<br>
>SKBR3_SR_sort_MAX.bam<br>


## 2. Subsampling
`samtools view --subsample` can perform subsampling of the data. First, the coverages of the "MAX" depth files must be calculated. 
```
#!/bin/bash

max_ONT=$(python get_genome_depth.py SKBR3_ONT_sort_MAX.bam)
max_PBCLR=$(python get_genome_depth.py SKBR3_PBCLR_sort_MAX.bam)
max_SR=$(python get_genome_depth.py SKBR3_SR_sort_MAX.bam)
```
Then calculate the `-s` or `--subsample` float, which is the fraction of target depth over MAX depth. In this analysis, I subsampled at depths at 5X to 30X with an interval of 5X. 
```
for target in {5,10,15,20,25,30}
do
ONT_FRAC=$(bc -l <<< $target/$max_ONT)
PBCLR_FRAC=$(bc -l <<< $target/$max_PBCLR)
SR_FRAC=$(bc -l <<< $target/$max_SR)
samtools view -s ${ONT_FRAC} -b SKBR3_ONT_sort_MAX.bam > SKBR3_ONT_sort_${d}X.bam
samtools view -s ${PBCLR_FRAC} -b SKBR3_PBCLR_sort_MAX.bam > SKBR3_PBCLR_sort_${d}X.bam
samtools view -s ${SR_FRAC} -b SKBR3_SR_sort_MAX.bam > SKBR3_SR_sort_${d}X.bam
done
```
There will be several files after this step. Example of resulting files with file name patterns:
>SKBR3_ONT_sort_5X.bam<br>
>SKBR3_ONT_sort_10X.bam<br>
>SKBR3_PBCLR_sort_30X.bam<br>
>SKBR3_SR_sort_10X.bam<br>

Lastly, index the files. 
```
for bamfile in $(ls *.bam)
do
samtools index $bamfile -@ 16
done
```

## 3. Calling structural variants (SVs)
Now, structural variant callers can be run on our long-read data. Some callers offer filtering options for output SVs based on minimum supporting reads. I suggest not using these options and filter the output raw VCFs yourself instead using `bcftools`. Then you won't need to rerun callers every time you decide to change your criteria. In my analysis, I organized my output directories following this format `SVcaller/sample/output` or `cuteSV/SKBR3_ONT_sort_5X/output`. <br><br>
`xx` here is from a method of batch submission. Including an example under cuteSV.
### cuteSV
I used the suggested cuteSV parameters for ONT and PacBio CLR data. Read about cuteSV [here](https://github.com/tjiangHIT/cuteSV). 
```
#!/bin/bash
Tech=$(echo "xx" | awk -F'[_]' '{print $2}')

if [[ $Tech == "ONT" ]]
then
cuteSV xx.bam GRCh38.p13.genome.fa xx.vcf cuteSV/xx/ \
--max_cluster_bias_INS 100 \
--diff_ratio_merging_INS 0.3 \
--max_cluster_bias_DEL 100 \
--diff_ratio_merging_DEL 0.3
elif [[ $Tech == "PBCLR" ]]
then
cuteSV xx.bam GRCh38.p13.genome.fa xx.vcf cuteSV/xx/ \
--max_cluster_bias_INS 100 \
--diff_ratio_merging_INS 0.3 \
--max_cluster_bias_DEL 200 \
--diff_ratio_merging_DEL 0.3
fi
```
<br>Batch submission:
```
q="'"
for file in *.bam ## careful here if LR and SR BAMs are not separated. Run LR callers on LR data.
do
name=$(basename ${file} .bam)
echo cp cuteSV.sh cuteSV_${file}.sh >> cuteSV_batch_submit.sh;
echo perl -pi -w -e ${q}s/xx/${name}/g${q} cuteSV_${file}.sh >> cuteSV_batch_submit.sh;
echo sbatch cuteSV_${file}.sh >> cuteSV_batch_submit.sh
mkdir cuteSV/${name}
done

# then run cuteSV_batch_submit.sh to the scheduler
```
### DeBreak
DeBreak includes a filter based on depth if given `--depth` option. Read about DeBreak [here](https://github.com/Maggi-Chen/DeBreak). An example using depth:
```
DepthX=$(echo "xx" | awk -F'[_]' '{print $4}')
Depth=$(STRING=$DepthX ; echo "${STRING//[!0-9]/}")

debreak --bam xx.bam \
-o DeBreak/xx/ -t 8 --rescue_large_ins --rescue_dup --poa --ref GRCh38.p13.genome.fa --depth $Depth
```
### PBSV
PBSV has two steps: `discover` and `call`. Use of the tandem repeats file is highly recommended. Read about PBSV [here](https://github.com/PacificBiosciences/pbsv).
```
pbsv discover --tandem-repeats human_GRCh38_no_alt_analysis_set.trf.bed xx.bam PBSV/xx/xx.svsig.gz
pbsv call -j 8 GRCh38.p13.genome.fa PBSV/xx/xx.svsig.gz PBSV/xx/xx.var.vcf
```
### Sniffles
Sniffles can also take advantage of a tandem repeats file. `minsvlen` is used here to just get SVs greater than or equal to 50 bp in length. Read about Sniffles [here](https://github.com/fritzsedlazeck/Sniffles).
```
sniffles --threads 4 -i xx.bam -v Sniffles/xx/xx.vcf \
--tandem-repeats human_GRCh38_no_alt_analysis_set.trf.bed \
--reference GRCh38.p13.genome.fa \
--minsvlen 50
```
### SVDSS
There are multiple steps to running SVDSS. Read about SVDSS [here](https://github.com/Parsoa/SVDSS).
```
# do this once
SVDSS index --fastq GRCh38.p13.genome.fa --index SVDSS/GRCh38.p13.genome.bwt --threads 16
```
```
# in a batch submission style
SVDSS smooth --bam xx.bam --workdir SVDSS/xx/ --reference GRCh38_p13/GRCh38.p13.genome.fa --threads 16

samtools index SVDSS/xx/smoothed.selective.bam SVDSS/xx/smoothed.selective.bam.bai -@ 16

SVDSS search --index SVDSS/GRCh38.p13.genome.bwt --bam SVDSS/xx/smoothed.selective.bam --workdir SVDSS/xx/ --assemble --threads 16

# For the next step, `--batches` is how many `.sfs` files were output from previous steps

dir=SVDSS/xx/
N=$(ls -dq $dir*sfs* | wc -l)

SVDSS call --reference GRCh38.p13.genome.fa --bam SVDSS/xx/smoothed.selective.bam --workdir SVDSS/xx/ --batches $N --threads 16
```
### SVIM
Read about SVIM [here](https://github.com/eldariont/svim).
```
svim alignment SVIM/xx/ xx.bam GRCh38.p13.genome.fa
```
