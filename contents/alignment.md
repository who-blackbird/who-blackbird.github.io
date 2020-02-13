# 4. Alignment

In this section we will cover:

- [Background](#background)
- [Data](#data)
- [Working Directory](#workingdirectory)
- [Reads Quality Control](#readsQC)
- [Alignment](#alignment)
- [Alignment Quality Control](#alignmentQC)
- [Visualisation](#visualisation)

You will learn to:

- Take the basecalled sequences in FASTQ format and align them to a genome of reference
- Perform quality control and visualise the alignment

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/alignment_schema.png" alt="alignment" class="inline"/>

## Background {#background}

**Coronaviruses** are enveloped positive-sense (+) single-stranded RNA viruses that infect a variety of mammalian and avian hosts. Among the several coronaviruses that are pathogenic to humans, most are associated with mild clinical symptoms, with the exception of severe acute respiratory syndrome (SARS) coronavirus (SARS-CoV), a Middle East respiratory syndrome (MERS) coronavirus (MERS-Cov), and the recently identified human-infecting coronavirus, provisionally named 2019 novel coronavirus (2019-nCoV, now Covid-19).

In this practical we will be using data from **Human Coronavirus 229E (HCoV-229E)**, one of the first coronavirus strains being described, with a genome size of ~27,300nt. In HCov-229E-infected cells, a total of seven major viral RNAs are produced. In it's 5'-terminal region, the genome RNA contains two large ORFs, 1a and 1b, that encode the viral replicase polyproteins 1a and 1b. mRNAs 2, 4, 5, 6 and 7 are used to produce the S protein, accessory protein 4, E protein, M protein and N protein respectively.

The scheme of genomic and subgenomic RNAs produced in HCov-229E-infected cells looks like this:

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/HCov-RNAs_CoV.png" alt="rnas" class="inline"/>

## Data {#data}

The data we will be using is from Human Coronavirus 229E (HCoV-229E), one of the first coronavirus strains being described.

In [this](http://www.ncbi.nlm.nih.gov/pubmed/?term=Viehweger+coronavirus) study, RNA sequencing was performed from Huh7 cells infected with serially passaged recombinant human coronaviruses, that we will refer to as:

- **WT RNA** (wild-type HCoV-229E)
- **SL2 RNA** (pool of recombinant virus, with the stem-loop structure (SL2) in the 5'UTR replaced with the equivalent SL2 element from SARS-CoV and SARS-BCoV)

Sequencing was performed with Oxford Nanopore MinION using the Direct RNA Sequencing protocol (SQK-RNA-1) and R9.4 chemistry. The source for the original FAST5 and FASTQ files can be found [here](https://osf.io/up7b4/).

In this practical we will take these FASTQ files and align them to the consensus reference genome of HCov-229E ([NC_002645.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_002645.1?report=genbank)). Then, we will assess the quality of the alignment and visualise it.

First of all, let's explore the FASTQ files!

A FASTQ file normally uses four lines per sequence:

1. Begins with a '@' and is followed by a sequence identifier
2. Is the raw sequence letters
3. Begins with a '+' character
4. Encodes the quality values for the sequence in Line 2

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/fastq.png" alt="img_1" class="inline"/>

To visualise the WT FASTQ file, you can type:

```
less -S ~/Course_Materials/data/alignment/fastq/WT_CoV.fastq.gz
```

and the SL2:

```
less -S ~/Course_Materials/data/alignment/fastq/SL2_CoV.fastq.gz
```

## Working Directory {#workingdirectory}

Now, we will set up the working directory where we will do the analysis. Open your terminal, go to

```
cd ~/Course_Materials/wd/day1
```

and create the following directory structure:

```
mkdir stats
mkdir alignment
```

Define now the following variables for convenience:

```
WT_fastq=~/Course_Materials/data/alignment/fastq/WT_CoV.fastq.gz
SL2_fastq=~/Course_Materials/data/alignment/fastq/SL2_CoV.fastq.gz
HCoV_ref=~/Course_Materials/data/alignment/reference_genome/HCov-229E.fasta
```

## Reads Quality Control {#readsQC}

There are many approaches to assess the quality of the reads. Here we will use [**NanoStat**](https://github.com/wdecoster/nanostat). This calculates various statistics from a long read sequencing dataset in FASTQ, BAM or albacore sequencing summary format.

For FASTQ files, it provides information for the number of reads, the read length and quality distributions (including mean and median) and the read length [N50](http://www.metagenomics.wiki/pdf/definition/assembly/n50).

From your wd, run:

```
NanoStat --fastq $WT_fastq > stats/WT_fastq_nanostats.txt
NanoStat --fastq $SL2_fastq > stats/SL2_fastq_nanostats.txt
```

Obtain for both WT and SL2 samples the:

- Mean and median read length
- Read length N50
- How many reads have a quality score (Q) higher than 15? For which sample?
- Longest read

_Hint:_ To visualize the file you can use the command `cat`, or `grep` if you want to search for lines containing a pattern, eg. `grep "Mean read length" stats/*fastq_nanostats.txt`

However, you may want to run your customised scripts, to answer your own questions... Here you have some ideas!

For example, to calculate how many reads we have in the FASTQ files you can also use `awk`:

```
zcat $WT_fastq | awk '{s++}END{print s/4}'
zcat $SL2_fastq | awk '{s++}END{print s/4}'
```

As well as to obtain the read length for each read:

```
zcat $WT_fastq  | awk '{if(NR%4==2) print length($1)}' > stats/WT_read_length.txt
zcat $SL2_fastq | awk '{if(NR%4==2) print length($1)}' > stats/SL2_read_length.txt
```

Then you can look at the read length distribution. For that, you can start R from the command-line:

```
R
```

and type the following:

```
#Load your libraries
library(ggplot2)
library(scales)

#Read your data
WTreadLength <- read.table("stats/WT_read_length.txt", header=FALSE, col.names = "length")
SL2readLength <- read.table("stats/SL2_read_length.txt", header=FALSE, col.names = "length")

#Add a column for the sample id
WTreadLength$sample <- "WT"
SL2readLength$sample <- "SL2"

#Merge the tables
readLength <- rbind(WTreadLength, SL2readLength)

#Reorder the groups
readLength$sample <- factor(readLength$sample, levels = c('WT', 'SL2'))

#Make the plot
p1 <- ggplot(data=readLength, aes(length, fill = sample)) +
    geom_histogram(position = "identity", alpha = 0.6, binwidth = 1) +
    scale_x_sqrt(breaks = trans_breaks("sqrt", function(x) x ^ 2)(c(1, max(readLength$length)))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Read length")
p1
```

You should see a plot like this one:

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/HCov_read_length.png" alt="readlength" class="inline"/>

For quitting R, just type:

```
quit()
```

## Alignment {#alignment}

Great! Now we will align these FASTQ files to the genome of reference - but first, a quick introduction about the format of the files you are going to handle.

The standard format for aligned sequence data is [SAM](http://samtools.github.io/hts-specs/SAMv1.pdf) (Sequence Alignment Map).

SAM files have a header that contains information on alignment and contigs used, and the aligned reads:

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/sam.jpg" alt="img_2" class="inline"/>

But because SAM files can be large, they are usually stored in the compressed version of them, [BAM](http://samtools.github.io/hts-specs/SAMv1.pdf) files.

Multiple algorithms have been developed to align long reads to a genome of reference. Some examples are:

- minimap2: [http://github.com/lh3/minimap2](http://github.com/lh3/minimap2)
- NGMLR: [http://github.com/philres/ngmlr](http://github.com/philres/ngmlr)
- LAST: [http://last.cbrc.jp](http://last.cbrc.jp)
- Graphmap: [http://github.com/isovic/graphmap](http://github.com/isovic/graphmap)
- bwa mem -x l ont2d: [http://github.com/lh3/bwa](http://github.com/lh3/bwa)

Here we will use **minimap2** to map the reads to the genome of reference. Then we will convert the SAM output to BAM format and sort it by mapping coordinate. For that we will use the following minimap2 options:

```
-a          Output in the SAM format (PAF by default)
-x splice   Long-read spliced alignment
-u f        How to find GT-AG. f:transcript strand
-k 14       K-mer size of 14 (as recommended)
```

**Note: The authors of the paper used the options ```-x map-ont``` (Nanopore vs reference mapping), ```--MD``` (Output the MD tag), and ```-u n``` (Don't match GT-AG) instead of ```-x splice -u f```. However, we updated them, since these are the current recommended optinos for Nanopore Direct RNA-seq**

We will start with the WT sample.

```
##Align to the ref using minimap
minimap2 -a -x splice -u f -k 14 $HCoV_ref $WT_fastq > alignment/WT_CoV.sam

##Sort and output as BAM
samtools sort alignment/WT_CoV.sam -O BAM -o alignment/WT_CoV.sort.bam
```

_Alternatively_, you can run these two steps using only one command line. Below is the alignment for the SL2 sample:

```
minimap2 -a -x splice -u f -k 14 $HCoV_ref $SL2_fastq | samtools sort - -O BAM -o alignment/SL2_CoV.sort.bam
```

Finally we will index the BAM files to run samtools subtools later.

```
samtools index alignment/SL2_CoV.sort.bam
samtools index alignment/WT_CoV.sort.bam
```

To visualise a BAM file:

```
samtools view -h alignment/SL2_CoV.sort.bam | less -S
```

## Alignment Quality Control {#alignmentQC}

As a first QC, we can run **NanoStat** on the BAM files:

```
NanoStat --bam alignment/WT_CoV.sort.bam > stats/WT_bam_nanostats.txt
NanoStat --bam alignment/SL2_CoV.sort.bam > stats/SL2_bam_nanostats.txt
```

Compare the:

- Mean and median aligned read length
- Number of reads aligned
- Read length N50
- Quality cutoffs
- Top highest quality and longest reads aligned

Do these results match with the values obtained from the fastq files?

Another option is to run **samtools stats**:

```
samtools stats alignment/WT_CoV.sort.bam > stats/WT_samtools_stats.txt
samtools stats alignment/SL2_CoV.sort.bam > stats/SL2_samtools_stats.txt
```

Explore the results obtained from samtools stats, and compare them with the NanoStat output.

- Which one provides more information?
- Compare the number of mapped vs unmapped reads. Why do you think there are so many unmapped ones?

Additional information for the `samtools stats` output can be found [here](https://www.htslib.org/doc/samtools-stats.html).

To obtain the coverage per base, we can use **samtools depth**:

```
samtools depth alignment/WT_CoV.sort.bam > stats/WT_samtools_depth.txt
samtools depth alignment/SL2_CoV.sort.bam > stats/SL2_samtools_depth.txt
```

And look at the coverage of both samples across the entire reference genome. For that, you can start R from the command-line:

```
R
```

and then, type the following:

```
#Load your libraries
library(ggplot2)

#Read your data
WTcov <- read.table("stats/WT_samtools_depth.txt", header=FALSE, col.names = c("chrom", "pos", "cov"))
SL2cov <- read.table("stats/SL2_samtools_depth.txt", header=FALSE, col.names = c("chrom", "pos", "cov"))


#Add a column for the sample id
WTcov$sample <- "WT"
SL2cov$sample <- "SL2"

#Merge the tables
coverage <- rbind(WTcov, SL2cov)

#Make the plot
p2 <- ggplot(data=coverage, aes(x = pos, y = cov, colour = sample)) +
    geom_line() +
    scale_y_log10() +
    scale_x_continuous(breaks = seq(0, max(coverage$pos), 1000)) +
    scale_colour_manual(values = c("#FF7F2D", "#0077B1")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Coverage", x = "Genome position", y = "Coverage")

p2
```

You should see a plot like this one.

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/HCoV-coverage.png" alt="cov" class="inline"/>

Which is a reproduction of the **Figure 2** from the original [paper](http://genome.cshlp.org/cgi/pmidlookup?view=long&pmid=31439691).

Additionally, you can also look at the coverage distribution in R:

```
library(dplyr)
library(purrr)
percent_cov <- coverage %>%
    group_by(sample) %>%
    mutate(percent = (map_int(cov, ~ sum(cov >= .x)))/27317) %>%
    select(sample, cov, percent) %>%
    unique

#Make the plot
p3 <- ggplot(data = percent_cov, aes(x = cov, y = percent*100, colour = sample)) +
    geom_line() +
    scale_x_continuous(breaks=seq(0,max(coverage$cov), 1000)) +
    labs(title = "Coverage distribution", x = "Coverage", y = "Percentage of bases")
p3
```

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/HCoV-coverage_distribution.png" alt="cov" class="inline"/>

You can also add a fancy vertical line to the previous plot intercepting with a minimum coverage of (let's say) 1000x:

```
p3 + geom_vline(xintercept = 1000, colour = "red")
```

## Visualisation {#visualisation}

To inspect the alignment, we will use [**Integrative Genomics Viewer**](https://software.broadinstitute.org/software/igv/).

Open IGV and load the HCov-229E.fasta reference genome by selecting Genomes>Load from File or Genomes>Load from URL. The new genome will be added to the drop-down menu, and also loaded and displayed in the IGV window.

Then, load your BAM file located in alignment/SL2_CoV.sort.bam

Explore!
