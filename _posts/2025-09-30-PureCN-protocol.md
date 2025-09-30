---
layout: post
title: "PureCN protocol"
date: 2025-09-30
---


##	This is the note how to run PureCN for mouse whole genome data

##	Background

Dr. Arun Pandiri has a cancer study involvoing a chemical call alpha pinene, I have helped with majority of analysis including

1. The sequenicing run was done at Sanger Institute, there, my analysis starts from the .bam file
2. I have processed the bam files, but removing low quality mapping, recalibrate the mapping quality
3. Since he does not have matched control, I constructed panel of normal for both male and female normal samples
4. Call the varian with Mutect2
5. Signature analysis with sigProfiler
6. Copy number analysis with cnvkit
7. Cancer driver gene analysis with dndscv
8. Strand bias -- not done

In order to evaluate and report cancer driver gene result we would like to apply the following four criteria

1. minimum two variants per gene
2. dndscv p-value < 0.05
3. must have impact on protein level (altered amino acid)
4. CCF > 30%

In order to calculate CCF, we need the purity and ploridy informaiton. From my literature search
I decide to use PureCN for this for the following benefit:

 
##	Available data

cnvkit output, segmentation file:

/Users/jyli/SciomeProjects/CCFstuff/AlphaPinene/cnvkitResults/FemaleLiverTumor/batchRun/
<table>
  <tr>
    <th>Owner</th>
    <th>Size</th>
    <th>File Name</th>
  </tr>
  <tr>
    <td>jyli</td>
    <td>48315135</td>
    <td>MD6818a_FemaleLiver_Tumor-diagram.pdf</td>
  </tr>
  <tr>
    <td>jyli</td>
    <td>75221</td>
    <td>MD6818a_FemaleLiver_Tumor-scatter.png</td>
  </tr>
  <tr>
    <td>jyli</td>
    <td>31</td>
    <td>MD6818a_FemaleLiver_Tumor.antitargetcoverage.cnn</td>
  </tr>
  <tr>
    <td>jyli</td>
    <td>4192233</td>
    <td>MD6818a_FemaleLiver_Tumor.bintest.cns</td>
  </tr>
  <tr>
    <td>jyli</td>
    <td>4019</td>
    <td>MD6818a_FemaleLiver_Tumor.call.cns</td>
  </tr>
  <tr>
    <td>jyli</td>
    <td>122013639</td>
    <td>MD6818a_FemaleLiver_Tumor.cnr</td>
  </tr>
  <tr>
    <td>jyli</td>
    <td>11251</td>
    <td>MD6818a_FemaleLiver_Tumor.cns</td>
  </tr>
  <tr>
    <td>jyli</td>
    <td>98136717</td>
    <td>MD6818a_FemaleLiver_Tumor.targetcoverage.cnn</td>
  </tr>

</table>


To get segmentation file that PureCN needs, I need to do the following two steps:

conda activate py38_cnvkit

cd /Users/jyli/SciomeProjects/CCFstuff/AlphaPinene/cnvkitResults/
cnvkit.py export seg  FemaleLiverTumor/batchRun/MD6818a_FemaleLiver_Tumor.cns -o MD6818a_FemaleLiver_Tumor.seg



##      Prepare binned genome and exons files

###	For genome, I have mm10.fa already

samtools faidx mm10.fa
cut -f1-2 mm10.fa.fai > mm10.genome

bedtools makewindows -g mm10.genome -w 100000 > mm10.100kb.bed

	If I ever want to create coverage, I can use the following commands
	But, it is unrealistics since the bam file is ~ 130 GB
	My laptop can't handle such a big file

	bedtools coverage -a mm10.100kb.bed -b my.bam > coverage_100kb.bed


###      To create the interval file
      Create interval file with PureCN with IntervalFile.R

Rscript $(Rscript -e "cat(system.file('extdata', 'IntervalFile.R', package='PureCN'))")  --infile mm10_10kb_bins_clean.bed   --fasta mm10.fa   --outfile  mm10_10kb_bins_cleaned.txt


<blockquote>
For exome, I want to create 100 bp window
It sounds straightforward, but I encountered 
Difficulty
</blockquote>

<blockquote>
	Failed!!!
	this crashes on mac os
curl -O http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/refGene.exonAll.bed.gz
gunzip refGene.exonAll.bed.gz
mv refGene.exonAll.bed mm10_exons.bed
</blockquote>




	I have to download Mus_musculus.GRCm38.102.gtf.gz
	And, it works

Download Ensembl GTF (replace version if needed)

curl -O ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz
gunzip Mus_musculus.GRCm38.102.gtf.gz


Extract only exon entries as BED

awk '$3 == "exon" {print $1"\t"$4-1"\t"$5}' Mus_musculus.GRCm38.102.gtf > mm10_exons.bed

bedtools sort -i mm10_exons.bed | bedtools merge -i - > mm10_exons_merged.bed
bedtools makewindows -b mm10_exons_merged.bed -w 100 > mm10_exon_bins_100bp.bed


Since I do NOT have the bam files, skip this
bedtools coverage -a mm10_exon_bins_100bp.bed -b my.bam > exon_bin_coverage_100bp.txt


Create interval file with PureCN with IntervalFile.R

###	For exome

Rscript $(Rscript -e "cat(system.file('extdata', 'IntervalFile.R', package='PureCN'))")  --infile mm10_exon_bins_100bp.bed   --fasta mm10.fa   --outfile mm10_exon_bins_100bp.txt

It calculates the gc_bias also, but no mappability
<pre>
Target	gc_bias	mappability	reptiming	Gene	on_target
chr1:3073253-3073613	0.279778393351801	NA	NA	.	TRUE
chr1:3073614-3073975	0.378453038674033	NA	NA	.	TRUE
chr1:3073976-3074337	0.364640883977901	NA	NA	.	TRUE
</pre>
Now, let's try to get exon with gene annotation

Extract exons with gene name
awk '$3 == "exon" {print $1"\t"$4-1"\t"$5"\t"$20}' Mus_musculus.GRCm38.102.gtf   | sed 's/"//g; s/;//g' > mm10_exons_with_gene.bed

bedtools makewindows -b mm10_exons_with_gene.bed -w 100 -i src > mm10_exon_bins_100bp_gene.bed
cat  mm10_exon_bins_100bp_gene.bed | grep -v "^JH" | grep -v "^GL" > temp  
mv temp mm10_exon_bins_100bp_gene.bed

Sort by gene then merge

bedtools sort -i mm10_exons_with_gene.bed   | bedtools groupby -g 4 -c 1,2,3 -o collapse > gene_exons_grouped.txt


Rscript $(Rscript -e "cat(system.file('extdata', 'IntervalFile.R', package='PureCN'))")  --infile mm10_exon_bins_100bp_gene.bed   --fasta mm10.fa   --outfile mm10_exon_bins_100bp_gene.txt


Unfortunately, I still miss gene information

cat mm10_exon_bins_100bp_gene.txt | grep -v ^@ | head 
<pre>
Target	gc_bias	mappability	reptiming	Gene	on_target
chr1:3073253-3073608	0.280898876404494	NA	NA	.	TRUE
chr1:3073609-3073965	0.380952380952381	NA	NA	.	TRUE
chr1:3073966-3074322	0.364145658263305	NA	NA	.	TRUE
chr1:3102016-3102125	0.363636363636364	NA	NA	.	TRUE
chr1:3205901-3206254	0.398305084745763	NA	NA	.	TRUE
</pre>
	Now, let' run it with transcriptome database
	Still NOT working, will NOT spend time on this
	Maybe until later
	The following line of command fails

Rscript $(Rscript -e "cat(system.file('extdata', 'IntervalFile.R', package='PureCN'))")  --infile mm10_exon_bins_100bp.bed   --fasta mm10.fa   --outfile mm10_exon_bins_100bp_withGeneSymbols.txt --txdb TxDb.Hsapiens.UCSC.hg19.knownGene


      Now, let's run PureCN



