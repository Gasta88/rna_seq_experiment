# rna_seq_experiment
This pages is for documenting the workflow that I've developed to compare [Kallisto](https://pachterlab.github.io/kallisto/about.html), [Salmon](https://combine-lab.github.io/salmon/), [TopHat](https://ccb.jhu.edu/software/tophat/index.shtml) and [STAR aligner](https://github.com/alexdobin/STAR).

This was the first project that I've made for the Data Analysis Group as a part-time member at the University of Dundee.

Raw data
========
The raw data comes from a [HT RNA-sequencing experiment](https://doi.org/10.1093/bioinformatics/btv425) run on S.cerevisiae with 7 technical replicates and 48 biological ones (336 in total) under two conditions (WT and Snf2).

A sample mapping is provided [here](https://figshare.com/articles/Metadata_for_a_highly_replicated_two_condition_yeast_RNAseq_experiment_/1416210).

Alignment
=========

Kallisto, Salmon and STAR jobs were sent to the School of Life Sciences cluster. All FASTQ files were stored on it in a filesystem that resemble the schema below:

```
~/dyeswap_data
	|---------/Snf2
		  |----/Snf2-1
	          |----/Snf2-2
	          |----/Snf2-3
	          |----/Snf2-4
	          |----/Snf2-5
	          |----/Snf2-6
	|---------/WT
		  |----/WT1
		  |----/WT2
		  |----/WT3
		  |----/WT4
		  |----/WT5
                  |----/WT6
```

Each tool need an index (Kallisto and Salmon) or a database (STAR) to correctly map the reads. The latest genome build was retrieved from Ensembl and used to create indexes for Kallisto and Salmon. 

The STAR database was already built and made available on the cluster at the School of Life Science.

Different Bash scripts were developed to run the applications on the cluster.

For Kallisto:

```bash
#!/bin/bash -l

for i in {1..6}
do
	mkdir out/output_Snf2-${i}
	fastqfiles=($(ls /homes/fgastaldello/dyeswap_data/Snf2/Snf2-${i}|egrep '\.fastq\.gz$'))
	module load kallisto&&kallisto quant -i SCtranscripts.idx -o out/output_Snf2-${i} -b 100 /homes/fgastaldello/dyeswap_data/Snf2/Snf2-${i}/${fastqfiles[0]} /homes/fgastaldello/dyeswap_data/Snf2/Snf2-${i}/${fastqfiles[1]}&&module unload kallisto
done

for i in {1..6}
do
        mkdir out/output_WT${i}
        fastqfiles=($(ls /homes/fgastaldello/dyeswap_data/WT/WT${i}|egrep '\.fastq\.gz$'))
        module load kallisto&&kallisto quant -i SCtranscripts.idx -o out/output_WT${i} -b 100 /homes/fgastaldello/dyeswap_data/WT/WT${i}/${fastqfiles[0]} /homes/fgastaldello/dyeswap_data/WT/WT${i}/${fastqfiles[1]}&&module unload kallisto
done
```

For Salmon:

```bash
#!/bin/bash -l

for i in {1..6}
do
	mkdir out/output_Snf2-${i}
	fastqfiles=($(ls /homes/fgastaldello/dyeswap_data/Snf2/Snf2-${i}|egrep '\.fastq\.gz$'))
	module load salmon&&salmon quant -i SCtranscript_index -l A -1 /homes/fgastaldello/dyeswap_data/Snf2/Snf2-${i}/${fastqfiles[0]} -2 /homes/fgastaldello/dyeswap_data/Snf2/Snf2-${i}/${fastqfiles[1]} -o out/output_Snf2-${i}&&module unload salmon
done

for i in {1..6}
do
        mkdir out/output_WT${i}
        fastqfiles=($(ls /homes/fgastaldello/dyeswap_data/WT/WT${i}|egrep '\.fastq\.gz$'))
        module load salmon&&salmon quant -i SCtranscript_index -l A -1 /homes/fgastaldello/dyeswap_data/WT/WT${i}/${fastqfiles[0]} -2 /homes/fgastaldello/dyeswap_data/WT/WT${i}/${fastqfiles[1]} -o out/output_WT${i}&&module unload salmon
done
```

For STAR:

```bash
#!/bin/bash -l

mkdir /homes/fgastaldello/STAR_test/out
for i in {1..6}
do
	mkdir /homes/fgastaldello/STAR_test/out/output_Snf2-${i}
	mkdir /homes/fgastaldello/STAR_test/out/output_WT${i}
done
module load "STAR/2.5"
for i in {1..6}
do
	echo ${i}
        cd /homes/fgastaldello/STAR_test/out/output_Snf2-${i}
	fastqfiles=($(ls /homes/fgastaldello/dyeswap_data/Snf2/Snf2-${i}|egrep '\.fastq\.gz$'))
	STAR --quantMode GeneCounts --runThreadN 10 --genomeDir /homes/fgastaldello/STAR_test/S.Cer_db --genomeLoad NoSharedMemory --readFilesCommand gunzip -c --readFilesIn /homes/fgastaldello/dyeswap_data/Snf2/Snf2-${i}/${fastqfiles[0]} /homes/fgastaldello/dyeswap_data/Snf2/Snf2-${i}/${fastqfiles[1]} --outSAMmode Full --outSJfilterIntronMaxVsReadN 50 100 150 200 --outFilterType BySJout --outFilterMultimapNmax 2 --outFilterMismatchNmax 5 --outSAMunmapped Within
done
for i in {1..6}
do
        echo ${i}
	cd /homes/fgastaldello/STAR_test/out/output_WT${i}
	fastqfiles=($(ls /homes/fgastaldello/dyeswap_data/WT/WT${i}|egrep '\.fastq\.gz$'))
	STAR --quantMode GeneCounts --runThreadN 10 --genomeDir /homes/fgastaldello/STAR_test/S.Cer_db --genomeLoad NoSharedMemory --readFilesCommand gunzip -c --readFilesIn /homes/fgastaldello/dyeswap_data/WT/WT${i}/${fastqfiles[0]} /homes/fgastaldello/dyeswap_data/WT/WT${i}/${fastqfiles[1]} --outSAMmode Full --outSJfilterIntronMaxVsReadN 50 100 150 200 --outFilterType BySJout --outFilterMultimapNmax 2 --outFilterMismatchNmax 5 --outSAMunmapped Within
done
```

Data analysis
=============

A simple data analysis was performed in [R language](https://www.r-project.org) using SpikeNorm and [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html). The former is a private package developed within the Data Analysis Group. It is developed by Nick Schurch and Pieta Schofield and it is used 
for processing RNA-seq data enriched with ERRC synthetic spike-ins transcripts.

The workflow is available in the `RNA-seq.Rmd` file in this repository. 

Conclusions
===========

This simple process is not perfect, but it was my first attempt in getting my hand dirty with RNA-seq experiement. The only conclusions that I can give are from a personal point of view.

Since I've never done something like this before, I have to say that user-friendliness in Kallisto and Salmon is greater than in STAR. 

The original data was mapped with TopHat, since it was considered a gold standard at that time. Personally, I found Kallisto and Salmon to be quicker and easier to use, in both setup and running time. 
In comparison STAR aligner has been more difficult to operate and slower to compute the mapping.

