# The Smithlab DNA Methylation Data Analysis Pipeline (MethPipe)

**Authors: Qiang Song, Benjamin Decato, Michael Kessler, Fang Fang, Jenny Qu, Tyler Garvin, Meng Zhou, Andrew Smith**

# Table of contents
1. [Assumptions](#assumptions)
2. [Citation information](#Citation_information)
3. [Methylome construction](#Methylome_construction)
	1. [Finding a Data Set](#Finding_a_Data_Set)
	2. [Mapping reads](#Mapping_reads)
	3. [Merging libraries and removing duplicates](#Merging_libraries_and_removing_duplicates)
4. []

<div style="text-align: justify"> 

The methpipe software package is a comprehensive pipeline and set of tools for analyzing whole genome bisulfite sequencing data (WGBS). This manual explains the stages in our pipeline, how to use the analysis tools, and how to modify the pipeline for your specific context.
</div>

## Assumptions <a name="assumptions"></a>

Our pipeline was designed to run in a cluster computing context, with many processing nodes available, and a jobsubmission system like PBS, SGE or SLURM, but it is also possible to analyze methylomes from most genomes(including human and mouse) in a local machine with at least 16 GB of RAM. Typically the data we deal withamounts to a minimum of 100GB for a mammalian methylome at 10x coverage. Intermediate files may cause thisamount to more than double during execution of the pipeline, and likely at the end of the pipeline the total size of fileswill amount to almost double the size of the raw data.Users are assumed to be quite familiar with UNIX/Linux and related concepts (e.g. building software from source,using the command line, shell environment variables, etc.).It is also critical that users are familiar with bisulfite sequencing experiments, especially the bisulfite conversionreaction, and how this affects what we observe in the sequenced reads. This is especially important if paired-end sequencingis used. If you do not understand these concepts, you will likely run into major problems trying to customizeour pipeline.

## Citation information <a name="Citation_information"></a>

If you find our programs helpful, please cite us. The majority of the tools described in this manual were introducedas part of the main methpipe manuscript, listed first below. However, several improvements and additions have beenmade over the years. Below is the main citation along with the papers describing those improvements or additions.
* methpipe: Qiang Song, Benjamin Decato, Elizabeth E Hong, Meng Zhou, Fang Fang, Jianghan Qu, Tyler Garvin, Michael Kessler, Jun Zhou, and Andrew D Smith. A reference methylome database and analysis pipeline to facilitate integrative and comparative epigenomics. PloS one, 8(12):e81148, 2013
* amrfinder: Fang Fang, Emily Hodges, Antoine Molaro, Matthew Dean, Gregory J Hannon, and Andrew D Smith. Genomic landscape of human allele-specific dna methylation. Proceedings of the National Academy of Sciences, 109(19):7332–7337, 2012* MLML: Jianghan Qu, Meng Zhou, Qiang Song, Elizabeth E Hong, and Andrew D Smith. Mlml: Consistent simultaneous estimates of dna methylation and hydroxymethylation. Bioinformatics, 29(20):2645–2646, 2013* Radmeth: Egor Dolzhenko and Andrew D Smith. Using beta-binomial regression for high-precision differential methylation analysis in multifactor whole-genome bisulfite sequencing experiments. BMC bioinformatics,15(1):1–8, 2014* PMD: Benjamin E Decato, Jianghan Qu, Xiaojing Ji, Elvin Wagenblast, Simon RV Knott, Gregory J Hannon, and Andrew D Smith. Characterization of universal features of partially methylated domains across tissues andspecies. Epigenetics & Chromatin, 2020

## Methylome construction <a name="Methylome_construction"></a>

### Finding a Data Det <a name="Finding_a_Data_Set"></a>

The key point of analyzing the data process is getting access to them. If you are a new MethPipe user at the very beginning of the data analysis path, we encourage you to go through the pipeline with a chosen at this moment dataset. One of the most popular databases storing gene expression datasets is [GEO](https://www.ncbi.nlm.nih.gov/geo/) (an abbreviation standing for Gene Expression Omnibus). Advanced searching on [this website](https://www.ncbi.nlm.nih.gov/gds/) provides you immediate access to numerous datasets. Start with your search with the above keywords:

* whole-genome bisulfite sequencing
* methylation profiling

**Quality Control**
Performing quality control on raw sequence data coming from high throughput sequencing pipelines is a required step in this type of analysis.  Checking, among others, GC content and sequence length distribution should become your new routine. We encourage you to use a faster version of FastQC software called Falco (FastQC Alternative Code). Falco is available on GitHub [here](https://github.com/smithlabcode/falco). 

### Mapping reads <a name="Mapping_reads"></a>
During bisulfite treatment, unmethylated cytosines in the original DNA sequences are converted to uracils, which are then incorporated as thymines (T) during PCR amplification. These PCR products are referred to as T-rich sequences as a result of their high thymine constitution. With paired-end sequencing experiments, the compliments of these T-rich sequences are also sequenced. These complimentary sequences have high adenine (A) constitution (A is the complimentary base pair of T), and are referred to as A-rich sequences. Mapping consists of finding sequence similarity, based on context specific criteria, between these short sequences, or reads, and an orthologous reference genome.When mapping T-rich reads to the reference genome, either a cytosine (C) or a thymine (T) in a read is considered a valid match for a cytosine in the reference genome. For A-rich reads, an adenine or a guanine is considered a valid match for a guanine in the reference genome. The mapping of reads to the reference genome by abismal is described below. If you choose to map reads with a different tool, make sure that your post-mapping files are appropriately formatted for the next components of the methpipe pipeline (necessary file formats for each step are covered in the corresponding sections). The default behavior of abismal is to assume that reads are T-rich and map accordingly, butdifferent sequencing protocols that generate A-rich and T-rich reads in different combinations are equally accepted. abismal is available on github [here](http://github.com/smithlabcode/abismal).

**Input and output file formats:** We assume that the original data is a set of sequenced read files, typically as producedby Illumina sequencing. These are FASTQ format files, and can be quite large. After the reads are mapped,these files are not used by our pipeline. The reference genome should be a single FASTA file that was previouslyindexed using the abismalidx tool. The abismal program requires an indexed FASTA reference genome and theinput FASTQ files(s), after which it generates a Sequence Alignment/Map(SAM) output indicating the coordinates ofmapped reads. Details of the SAM file format can be found at the SAM file format documentation. These SAM fileswill be the input files for the postprocessing quality control and analysis programs to follow, including bsrate andmethcounts.Abismal operates by preprocessing the reference genome into a large index, where k-mers of set length are usedas keys to a list of potential mapping locations for reads that begin with their sequence. To produce this index run thefollowing command:
<pre><code>abismalidx genome_folder_or_file index_file
</code></pre>

This index file is roughly 2.5 times larger than the input reference genome size. For the human genome, whosesize is 3 GB, the resulting index is approximately 7.5 GB.

**Decompressing and isolating paired-end reads:** Sometimes paired-end reads are stored in the same FASTQ file. Because we treat these paired ends differently, they must be separated into two files and both files must be passed asinputs to abismal. If your data is compressed as a Sequenced Read Archive, or SRA file, you can decompress and split pairedend reads into two files at the same time using fastq-dump, which is a program included in the sra-toolkitpackage, available for most unix systems. Below is an example of using fastq-dump to decompress and separateFASTQ data by end:<pre><code>$ fastq-dump --split-3 human_esc.sra</code></pre>

If you have a FASTQ file not compressed in SRA format, you can split paired ends into two separate files byrunning the following commands:
<pre><code>$ sed -ne ’1˜8{N;N;N;p}’ *.fastq > *_1.fastq
$ sed -ne ’4˜8{N;N;N;p}’ *.fastq > *_2.fastq
</code></pre>

**Sequencing adapters:** These are a problem in any sequencing experiment with short fragments relative to the lengths of reads. A robust method for removing adapters is available via the Babraham Institute’s Bioinformaticsgroup, called trim galore. This program takes into account adapter contamination of fewer than 10 nucleotides and can handle mismatches with the provided sequence. We strongly recommend reads are trimmed prior to mapping. Our recommended parameter choice for trim galore restricts trimming only to sequencing adapters, leaving all other bases unaltered. This can be attained by running it with the following command for single-end mapping<pre><code>$ trim_galore -q 0 --length 0 human_esc.fastq
</code></pre>

and the following command for paired-end<pre><code>$ trim_galore --paired -q 0 --length 0 human_esc_1.fastq human_esc_2.fastq
</code></pre>Note that the --paired flag is necessary to ensure the program does not interpret the two input files as independentand that the resulting FASTQ files still have corresponding read mates in the correct order.

**Single-end reads:** When working with data from a single-end sequencing experiment, you will have T-rich readsonly. Abismal expects T-rich reads as a default. Execute the following command to map all of your single-end readswith abismal:
<pre><code>$ abismal -i <genome index> -o output_SAM [options] input_fastq</code></pre>
To save files in BAM format, which significantly reduce disk space, simply redirect the abismal output to the samtools view program using the -b flag to compress to BAM and the -S flag to indicate that input is in SAM format.<pre><code>$ abismal -i genome_index -m output_STATS [options] input\_fastq |
samtools view -bS > output_BAM
</code></pre>

**Paired-end reads:** When working with data from a paired-end sequencing experiment, you will have T-rich and A-rich reads. T-rich reads are often kept in files labeled with an “ 1” and A-rich reads are often kept in files labeledwith an “ 2”. T-rich reads are sometimes referred to as 50 reads or mate 1 and A-rich reads are sometimes referred to 30 reads or mate 2. We assume that the T-rich file and the A-rich contain the same number of reads, and eachpair of mates occupy the same lines in their respective files. We will follow this convention throughout the manual and strongly suggest that you do the same. Run the following command to map two reads files from a paired-endsequencing experiment:<pre><code>$ abismal -i index -o output_SAM [options] input_fastq 1 input_fastq_2
</code></pre>

In brief, what happens internally in abismal is as follows. abismal finds candidate mapping locations for a T-rich mate with CG-wildcard mapping, and candidate mapping locations for the corresponding A-rich mate withAG-wildcard mapping. If two candidate mapping locations of the pair of mates are within certain distance in the same chromosome and strand and with correct orientation, the two mates are combined into a single read (after reversecomplement of the A-rich mate), referred to as a fragment. The overlapping region between the two mates, if any, is included once, and the gap region between them, if any, is filled with Ns. The parameters -l and -L to abismalindicate the minimum and maximum size of fragments to allow to be merged, respectively.Here the fragment size is the sum of the read lengths at both ends, plus whatever distance is between them. So this is the length of the original molecule that was sequenced, excluding the sequencing adapters. It is possible for a given read pair that the molecule was shorter than twice the length of the reads, in which case the ends of the mates will overlap, and so in the merged fragment will only be included once. Also, it is possible that the entire molecule was shorter than the length of even one of the mates, in which case the merged fragment will be shorter than either of the read ends. If the two mates cannot be merged because they are mapped to different chromosomes or differentstrand, or they are far away from each other, abismal will output each mate individually if its mapping position is unambiguous.abismal provides a statistical summary of its mapping job in the “mapstats” file, which includes the total numberand proportion of reads mapped, how many paired end mates were mapped together, and the distribution of fragment lengths computed by the matched pairs. The “mapstats” file name is usually the same as the SAM output with “.mapstats” appended to it, but custom file names can be provided using the -m flag.**Mapping reads in a large file:** Mapping reads may take a while, and mapping reads from WGBS takes even longer. It usually take quite a long time to map reads from a single large file with tens of millions of reads. If you have access to a cluster, one strategy is to launch multiple jobs, each working on a subset of reads simultaneously, and finally combine their output. abismal takes advantage of OpenMP to parallelize the process of mapping reads using theshared index. 
If each node can only access its local storage, dividing the set of reads to into k equal sized smaller reads files, and mapping these all simultaneously on multiple nodes, will make the mapping finish about k times faster. The UNIX split command is good for dividing the reads into smaller parts. The following BASH commands will take a directory named reads containing Illumina sequenced reads files, and split them into files containing at most 3Mreads:<pre><code>$ mkdir reads_split$ for i in reads/*.txt; do \split -a 3 -d -l 12000000 ${i} reads_split/$(basename $i); done
</code></pre>
Notice that the number of lines per split file is 12M, since we want 3M reads, and there are 4 lines per read. If you split the reads like this, you will need to “unsplit” them after the mapping is done. This can be done using thesamtools merge command.

**Formatting reads:** Before analyzing the output SAM generated by mappers, some formatting is required. The first formatting step is to merge paired-end mates into single-end entries. This is particularly important to quantify methylation, as fragments that overlap must count the overlapping bases only once and must be treated as originating from the same allele. These can be ensured by merging them into a single entry. In what follows we will use the human esc placeholder name to describe an example dataset to be analyzed. SAM files generated by abismal can be formatted using the format reads program as shown below:<pre><code> $ format_reads -o human_esc_f.sam -f abismal human_esc.bam</code></pre>
Here we added the “ f” suffix to indicate that the SAM file was generated after the format reads program was run. In addition to abismal described above, users may also wish to process raw reads using alternative mapping algorithms, including Bismark and BSMAP. These WGBS mappers similarly generate SAM files, so the program format reads can also be used to convert the output from those mappers to the standardized SAM format used in our pipeline. To convert BSMAP mapped read file in .bam format, run<pre><code>$ format_reads -o human_esc_f.sam -f bsmap human_esc.bam</code></pre>where the option -f specifies that the original mapper is BSMAP. To obtain a list of alternative mappers supported by our converter, run format reads without any options.

### Merging libraries and removing duplicates <a name="Merging_libraries_and_removing_duplicates"></a>
This is a sub paragraph, formatted in heading 3 style


### Sub paragraph <a name="subparagraph1"></a>
This is a sub paragraph, formatted in heading 3 style