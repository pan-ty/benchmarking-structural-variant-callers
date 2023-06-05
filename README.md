# benchmarking-structural-variant-callers
### Description
In this project, I benchmarked several long-read structural variant (SV) callers (cuteSV, DeBreak, PBSV, Sniffles, SVDSS, and SVIM) on a range of sequencing depths (5X, 10X, 15X, 20X, 25X, 30X, MAX) in order to determine the minimum depth allowing for the most sufficient yield of structural variant calls. With these benchmarking statistics, a cost-effective sequencing depth can be determined to generate new data for SV studies. For this analysis, I used data of the SK-BR-3 breast cancer cell line (SKBR3) from three sequencing platforms: Oxford Nanopore Technologies (ONT), PacBio Continuous Long Reads (PBCLR), and Illumina short-reads (SR).<br>

 ## 1. Alignment
 I used minimap2 to align long-reads (ONT and PBCLR) and bwa-mem2 to align the short-reads (SR) on GRCh38.p13.
```bash
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
```bash
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
```bash
max_ONT=$(python get_genome_depth.py SKBR3_ONT_sort_MAX.bam)
max_PBCLR=$(python get_genome_depth.py SKBR3_PBCLR_sort_MAX.bam)
max_SR=$(python get_genome_depth.py SKBR3_SR_sort_MAX.bam)
```
Then calculate the `-s` or `--subsample` float, which is the fraction of target depth over MAX depth. In this analysis, I subsampled at depths at 5X to 30X with an interval of 5X. 
```bash
for target in {5,10,15,20,25,30}
do
ONT_FRAC=$(bc -l <<< $target/$max_ONT)
PBCLR_FRAC=$(bc -l <<< $target/$max_PBCLR)
SR_FRAC=$(bc -l <<< $target/$max_SR)
samtools view -s ${ONT_FRAC} -b SKBR3_ONT_sort_MAX.bam > SKBR3_ONT_sort_${d}X.bam
samtools view -s ${PBCLR_FRAC} -b SKBR3_PBCLR_sort_MAX.bam > SKBR3_PBCLR_sort_${d}X.bam
samtools view -s ${SR_FRAC} -b SKBR3_SR_sort_MAX.bam > SKBR3_SR_sort_${d}X.bam
done
```bash
There will be several files after this step. Example of resulting files with file name patterns:
>SKBR3_ONT_sort_5X.bam<br>
>SKBR3_ONT_sort_10X.bam<br>
>SKBR3_PBCLR_sort_30X.bam<br>
>SKBR3_SR_sort_10X.bam<br>

Lastly, index the files. 
```bash
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
```bash
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
```bash
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
```bash
DepthX=$(echo "xx" | awk -F'[_]' '{print $4}')
Depth=$(STRING=$DepthX ; echo "${STRING//[!0-9]/}")

debreak --bam xx.bam \
-o DeBreak/xx/ -t 8 --rescue_large_ins --rescue_dup --poa --ref GRCh38.p13.genome.fa --depth $Depth
```
### PBSV
PBSV has two steps: `discover` and `call`. Use of the tandem repeats file is highly recommended. Read about PBSV [here](https://github.com/PacificBiosciences/pbsv).
```bash
pbsv discover --tandem-repeats human_GRCh38_no_alt_analysis_set.trf.bed xx.bam PBSV/xx/xx.svsig.gz
pbsv call -j 8 GRCh38.p13.genome.fa PBSV/xx/xx.svsig.gz PBSV/xx/xx.var.vcf
```
### Sniffles
Sniffles can also take advantage of a tandem repeats file. `minsvlen` is used here to just get SVs greater than or equal to 50 bp in length. Read about Sniffles [here](https://github.com/fritzsedlazeck/Sniffles).
```bash
sniffles --threads 4 -i xx.bam -v Sniffles/xx/xx.vcf \
--tandem-repeats human_GRCh38_no_alt_analysis_set.trf.bed \
--reference GRCh38.p13.genome.fa \
--minsvlen 50
```
### SVDSS
There are multiple steps to running SVDSS. Read about SVDSS [here](https://github.com/Parsoa/SVDSS).
```bash
# do this once
SVDSS index --fastq GRCh38.p13.genome.fa --index SVDSS/GRCh38.p13.genome.bwt --threads 16
```
```bash
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
```bash
svim alignment SVIM/xx/ xx.bam GRCh38.p13.genome.fa
```

### SV calling on short-read data
We are not interested in benchmarking short-read callers as our future studies involve only long-read data. It can still be useful in building a comprehensive groundtruth set. Read about [lumpy](https://github.com/arq5x/lumpy-sv) and [SvABA](https://github.com/walaj/svaba). 
```bash
# lumpy
lumpyexpress -B SKBR3_SR_sort_MAX.bam -o SKBR3_SR_sort_MAX.vcf

# SvABA
svaba run -t SKBR3_SR_sort_MAX.bam -p 8 -D dbsnp_indel.vcf -a SKBR3_SR_sort_MAX -G GRCh38.p13.genome.fa

```

## 4. Filtering VCFs
To filter the VCFs, it is important to first examine the information field given by each caller and determine the best filtering criteria. 
For my analysis, I filtered each depth according to the minimum supporting reads below. Based on my preliminary benchmarking statistics, the filtering used here can and should be improved upon. 

| Depth | Min Supp |
| ----- | -------- |
| 5X | 2 |
| 10X | 3 |
| 15X | 4 |
| 20X | 5 |
| 25X | 6 |
| 30X | 7 |
| 59X (SR MAX) | 12 | 
| 61X (ONT MAX) | 12 | 
| 112X (PBCLR MAX) | 15 | 

The following is an example of filtering all of the cuteSV VCFs. Here, the directory cuteSV contains several subdirectories, each named according to the file name pattern and contains the VCF output from the caller. `minsup_byname.py SKBR3_ONT_sort_5X` will output the minimum supporting reads for 5X, which is `2` based on the table above. Then `bcftools view -i` can filter the `INFO/RE` field using `$MinSup`, along with other necessary filters. The `filtervcf_based_on_length.py` filters based on SV length or the `SVLEN` field, only keeping calls with an absolute value greater than or equal to 50. This is important for accurate benchmarking since the ground truth set will only contain SVs >= 50 bp in length. 

### cuteSV
```bash
for name in cuteSV/*/
do
MinSup=$(python minsup_byname.py ${name})
cd cuteSV/${name}/
bcftools view -i 'FILTER == "PASS" && INFO/PRECISE == 1 && INFO/RE >= '$MinSup ${samplename}.vcf > temp_cuteSV_${name}.vcf
python filtervcf_based_on_length.py -i temp_cuteSV_${name}.vcf -o ~/filtered_vcf/cuteSV_${name}.vcf -l 50
rm temp_cuteSV_${name}.vcf
done
```
### other callers
The final location of the final filtered VCFs go into a folder titled `filtered_vcf` with the SV tool name added to the beginning of the file names. There is room for improvement on these filtering criteria.  
```bash
# DeBreak
bcftools view -i 'FILTER == "PASS" && INFO/PRECISE == 1 && INFO/SUPPREAD >= '"$MinSup"'' ${name}.vcf > temp_DeBreak_${name}.vcf
python filter_vcf_based_on_length.py -i temp_DeBreak_${name}.vcf -o ~/filtered_vcf/DeBreak_${name}.vcf -l 50

# PBSV
bcftools filter -i 'FILTER=="PASS" && FORMAT/AD[0:1]>='"${MinSup}"  ${name}.var.vcf > temp_PBSV_${name}.vcf
python filter_vcf_based_on_length.py -i temp_PBSV_${name}.vcf -o ~/filtered_vcf/PBSV_${name}.vcf -l 50

# Sniffles
bcftools view -i 'INFO/SUPPORT >= '"$MinSup"' && INFO/IMPRECISE != 1' ${name}.vcf > temp_Sniffles_${name}.vcf
python filter_vcf_based_on_length.py -i temp_Sniffles_${name}.vcf -o ~/filtered_vcf/Sniffles_${name}.vcf -l 50

# SVDSS
bcftools view -i 'FILTER == "PASS" && INFO/COV >= '"$MinSup" svs_poa.vcf > temp_SVDSS_${name}.vcf
python filter_vcf_based_on_length.py -i temp_SVDSS_${name}.vcf -o ~/filtered_vcf/SVDSS_${name}.vcf -l 50

# SVIM
bcftools filter -i 'QUAL >= '"$MinSup"'&& FILTER == "PASS" && INFO/SUPPORT >= 2' variants.vcf > temp_SVIM_${name}.vcf
python filter_vcf_based_on_length.py -i temp_SVIM_${name}.vcf -o ~/filtered_vcf/SVIM_${name}.vcf -l 50

# lumpy
# The VCFs should be reheaded with contigs first
bcftools reheader --fai GRCh38.p13.genome.fa.fai SKBR3_SR_sort_MAX.vcf -o rehead_lumpy_SKBR3_SR_sort_MAX.vcf
bcftools view -i 'INFO/PE >='"${MinSup}" rehead_lumpy_SKBR3_SR_sort_MAX.vcf > temp_lumpy_SKBR3_SR_sort_MAX.vcf
python filter_vcf_based_on_length.py -i temp_lumpy_SKBR3_SR_sort_MAX.vcf -o lumpy_SKBR3_SR_sort_MAX.vcf -l 50

# SvABA
# to convert BND to basic SV types on the already-filtered output .svaba.sv.vcf
python SVclassifier_SvABA.py -i SKBR3_SR_sort_MAX.svaba.sv.vcf -o ~/filtered_vcf/svaba_${name}.vcf
```

## 5. Preparing the ground truth set
To prepare the ground truth for benchmarking, I used the package `SURVIVOR` and the MAX depth VCFs from all callers, both long-read and short-read. Read about SURVIVOR [here](https://github.com/fritzsedlazeck/SURVIVOR/wiki)
First, I made sure the VCFs for `SURVIVOR` were sorted.
```bash
cd filtered_vcf
for maxfile in *MAX.vcf
do
bcftools sort $maxfile -o sorted_$maxfile
done
```
Then I ran `SURVIVOR merge` using intersections and unions as outlined below. In this analysis, intersections used options `50 2 0 0 0 50`, meaning that the maximum allowed distance between SVs was 50, the agreement of SV type and strand was disregarded, and only SVs greater or equal to 50 bp in length were compared. For unions, the `2` was changed to `0`. . 
```bash
# Intersections. Within ONT, PBCLR, SR
ls sorted*ONT*.vcf > ONT_sample_files
SURVIVOR merge ONT_sample_files 50 2 0 0 0 50 truth_set/ONT_LR_gt.vcf

ls sorted*PBCLR*.vcf > PBCLR_sample_files
SURVIVOR merge PBCLR_sample_files 50 2 0 0 0 50 truth_set/PBCLR_LR_gt.vcf

ls sorted*.vcf > SR_sample_files
SURVIVOR merge SR_sample_files 50 2 0 0 0 50 truth_set/SR_gt.vcf

# Union. ONT and PBCLR for LR ground truth
cd truth_set
ls *LR_gt.vcf > LR_sample_files
SURVIVOR merge LR_sample_files 50 0 0 0 0 50 union_LR_gt.vcf

# Union. LR and SR for final ground truth
ls union_LR_gt.vcf SR_gt.vcf > final_sample_files
SURVIVOR merge final_sample_files 50 0 0 0 0 50 final_gt.vcf
```

## 6. Benchmarking
`TRUVARI` was used to for benchmarking. This package requires all VCFs used in comparisons in the `.gzip` and `.gzip.tbi` formats. Read about TRUVARI [here](https://github.com/ACEnglish/truvari).
```bash
# gzip and index the filtered VCFs
for file in filtered_vcf/*.vcf
do
bcftools sort ${file} -Oz -o filtered_vcf/processed/${file}.gz
bcftools index -t filtered_vcf/processed/${file}.gz
done

# gzip and index ground truth VCF
bcftools sort truth_set/final_gt.vcf -Oz -o truth_set/final_gt.vcf.gz
bcftools index -t truth_set/final_gt.vcf.gz


# Run benchmarking
gt=truth_set/final_gt.vcf.gz
for sample in filtered_vcf/processed/*.vcf.gz
do
name=$(basename ${sample} .vcf.gz)
truvari bench -b ${gt} -c filtered_vcf/processed/${sample} -o benchmarking/${sample}
done
```

## 7. Results
The following are methods I used to extract benchmarking results, number of SV calls and types by caller, and ground truth information so I could more easily create plots in R. 
### Benchmarking statistics
`TRUVARI` outputs a `summary.json` file for each comparison with the ground truth (e.g., PBSV_SKBR3_ONT_sort_25X.vcf with final_gt.vcf). `get_benchmark.py` should create a table with the Caller, Tech (ONT/PBCLR), Depth, Precision, Recall, and F1 as columns, and corresponding values from each comparison in the rows. This only works if correct file name patterns and output directory were used. 
```bash
results=benchmarking/
# subdirectories created by TRUVARI within the benchmarking folder might look like cuteSV_SKBR3_ONT_sort_5X/ and DeBreak_SKBR3_PBCLR_sort_30X
python get_benchmark.py ${results}
```
### SV calls by depth and caller
Using `SURVIVOR stats` and bash version to get total SVs, DEL, DUP, INS, INV, TRA in a table:
```bash
# header for results file
echo -e "Caller\tTech\tDepth\tTot\tDEL\tDUP\tINS\tINV\tTRA" >> results/caller_results.txt

cd filtered_vcf

for file in *.vcf
do
name=$(basename ${file} .vcf)
DepthX=$(echo ${name} | awk -F'[_]' '{print $5}')
Tech=$(echo ${name} | awk -F'[_]' '{print $3}')
Caller=$(echo ${name} | awk -F'[_]' '{print $1}')

SURVIVOR stats ${name}.vcf vcf_summary -1 -1 -1 > results/stats/stats_${name}.txt

stat=$(sed -n '4,4p;4q' results/stats/stats_${name}.txt)
echo -e "${Caller}\t${Tech}\t${DepthX}\t${stat}" >> results/caller_results.txt
done
```
### Ground truth
A similar method using bash for a table of SVs contained by the ground truth VCF:
```bash
# header for results file
echo -e "Tot\tDEL\tDUP\tINS\tINV\tTRA" >> results/truth_results.txt

SURVIVOR stats truth_set/final_gt.vcf vcf_summary -1 -1 -1 > results/stats/stats_final_gt.txt
stat=$(sed -n '4,4p;4q' results/stats/stats_final_gt.txt)
echo -e "${stat}" >> results/truth_results.txt
```

After acquiring these tables, I created plots in R to visualize [the results]. 
