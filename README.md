# CRISPR-Cpf1_indels_analysis
A pipeline to analyze features of CRISPR/Cpf1 induced indels, including frequencies, sizes and positions

## Overview 
This pipeline analyzes frequencies, sizes and relative positions for CRISPR/Cpf1 induced small insertions/deletions (indels). It takes fastq files as inputs and outputs a list of sizes and positions for each detected indel and a metadata of overall frequencies of deletions and insertions separately and combined as well as total number of unique indels. 

## bwa_mapping.sh
**bwa_mapping.sh is the first script to be called.** <br />
**Inputs**: fastq files for targeted illumina sequencing of PCR amplified guide RNA targets and flanking regions <br />
**Outputs**: <br />
1. .sam file of mapping information 
2. .bam file of mapping information 
3. .csv file of sizes and positions 
4. .meta file (text) of overall deletions and insertions frequencies as well as total number of unique indels 

**Steps**: <br />
1. In the Harvard HPCC environment (orchestra), first load all necessary packages/softwares
2. Iterate over each fastq file and map to the indexed reference “genome” (guide RNA targets and flanking regions) using bwa mem (gap open penalty = 10; gap extension penalty = 1)
3. Convert the resulting .sam file from step 2 into sorted and indexed .bam file for igv visualization 
4. Call VariantsCalling.py to extract and process SIGAR string from the .sam files to calculate overall and individual indel frequency/size/position 

## VariantCalling.py
**VariantsCalling.py is called in bwa_mapping.sh.** <br />
**Input**: .sam file from bwa mapping <br />
**Output**: a .csv file for a list of indel size/position and a .metadata file for overall indel/insertion/deletion frequency and total number of unique indels for each sample <br />

## IndelAnalysis.R 
**IndelAnalysis.R is called after bwa_mapping.sh.** <br />
**Input**: .csv and .meta files from VariantsCalling.py <br />
**Output**:
1. A bar plot for overall indel frequency (deletions and insertions shown separately) 
2. A histogram of indel sizes 
3. A histogram of indel positions (relative position to 3’ guide sequence)
4. A histogram of deletions in the guide RNA sequence

A Rmarkdown file is generated to show the flow of the analysis for eGFP data. The input data and .Rmd files are included in this repository. The published .html file can be found here: http://rpubs.com/missingboy/350343
