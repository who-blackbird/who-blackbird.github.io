# 4. Alignment

In this section we will cover:

* [Data](#Data)
* [WorkingDirectory](#WorkingDirectory)
* [ReadsQC](#ReadsQC)
* [Alignment](#Alignment)
* [AlignemntQC](#AlignemntQC)
* [Visualisation](#Visualisation)

You will learn to:
- Take the basecalled sequences in FASTQ format and align them to a genome of reference
- Perform quality control and visualise the alignment

## 1. Data

The data we will be using is from Human Coronavirus 229E (HCoV-229E), one of the first coronavirus strains being described, with a genome size of ~27,300nt.

In [this](http://www.ncbi.nlm.nih.gov/pubmed/?term=Viehweger+coronavirus) study, RNA sequencing was performed from Huh7 cells infected with serially passaged recombinant human coronaviruses, that we will refer to as:

* WT RNA (wild-type HCoV-229E)
* SL2 RNA (pool of recombinant virus, with the stem-loop structrue (SL2) in the 5'UTR replaced with the equivalent SL2 element from SARS-CoV and SARS-BCoV)

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
less -S ~/Course_Materials/nanopore_practical/data/fastq/2017-09-05_coronavirus_WT.rna.fastq.gz
```

and the SL2:

```
less -S ~/Course_Materials/nanopore_practical/data/fastq/2017-09-29_coronavirus_SL2.rna.fastq.gz
```

## 2. Working directory

First we will set up the working directory where we will do the analysis. Open your terminal, go to 

```
cd ~/Course_Materials/nanopore_practical/wd
```

and create the following directory structure:

```
mkdir stats
mkdir alignment
mkdir variant_calling
mkdir annotation
```

## 2. Reads QC

We will now calculate how many reads we have in the FASTQ files:

```
awk '{s++}END{print s/4}' ~/Course_Materials/nanopore_practical/data/fastq/WT_CoV.fastq.gz
awk '{s++}END{print s/4}' ~/Course_Materials/nanopore_practical/data/fastq/SL2_CoV.fastq.gz
```

You can then use awk to obtain the read length for each read:

```
awk '{if(NR%4==2) print length($1)}' ~/Course_Materials/nanopore_practical/data/fastq/WT_CoV.fastq.gz > stats/WT_read_length.txt
awk '{if(NR%4==2) print length($1)}' ~/Course_Materials/nanopore_practical/data/fastq/SL2_CoV.fastq.gz > stats/SL2_read_length.txt
```

And look at the read length distribution. For that, you can start R from the command-line:

```
R
```

and then, type the following:

```
#Load your libraries
library(ggplot2)

#Read your data
WTreadLength <- read.table("stats/WT_read_length.txt", header=FALSE, col.names = "length")
SL2readLength <- read.table("stats/SL2_read_length.txt", header=FALSE, col.names = "length")

#Add a column for the sample id
WTreadLength$sample <- "WT"
SL2readLength$sample <- "SL2"

#Merge the tables
readLength <- rbind(WTreadLength, SL2readLength)

#Make the plot
ggplot(data=readLength, aes(length, fill = sample)) + 
    geom_histogram(color="#e9ecef", position = 'dodge', binwidth = 10) +
    scale_fill_manual(values=c("#69b3a2", "#404080"))
```

For quitting R, just type:

```
quit()
```

## Alignment

Great! Now we will align these FASTQ files to the genome of reference - but first, a quick introduction about the format of the files you are going to handle.

The standard format for aligned sequence data is [SAM](http://samtools.github.io/hts-specs/SAMv1.pdf) (Sequence Alignment Map).

SAM files have a header that contains information on alignment and contigs used, and the aligned reads:

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/sam.jpg" alt="img_2" class="inline"/>

But because SAM files can be large, they are usually stored in the compressed version of them, [BAM](http://samtools.github.io/hts-specs/SAMv1.pdf) files.

Multiple algorithms have been developed to align long reads to a genome of reference. Some examples are:
-	minimap2: [http://github.com/lh3/minimap2](http://github.com/lh3/minimap2)
-	NGMLR: [http://github.com/philres/ngmlr](http://github.com/philres/ngmlr)
-	LAST: [http://last.cbrc.jp](http://last.cbrc.jp)
-	Graphmap: [http://github.com/isovic/graphmap](http://github.com/isovic/graphmap)
-	bwa mem -x l ont2d: [http://github.com/lh3/bwa](http://github.com/lh3/bwa)

Here we will use **minimap2** to map the reads to the genome of reference. Then we will convert the SAM output to BAM format and sort it by mapping coordinate.

```
##Align to the ref using minimap
minimap2 -x map-ont --MD -u n -k 14 -a ~/Course_Materials/nanopore_practical/data/reference_genome/HCov-229E.fasta ~/Course_Materials/nanopore_practical/data/fastq/WT_CoV.fastq.gz > alignment/SL2_CoV.sam

##Sort and output as BAM
samtools sort alignment/SL2_CoV.sam -O BAM -o alignment/SL2_CoV.sort.bam
```

***Alternatively***, you can run these two steps using only one command line:

```
minimap2 -x map-ont --MD -u n -k 14 -a ~/Course_Materials/nanopore_practical/data/reference_genome/HCov-229E.fasta ~/Course_Materials/nanopore_practical/data/fastq/WT_CoV.fastq.gz | samtools sort - -O BAM -o alignment/SL2_CoV.sort.bam
```

Finally we will index the BAM file to run samtools subtools later.

```
samtools index alignment/SL2_CoV.sort.bam
```

To visualise the BAM file:

```
samtools view alignment/SL2_CoV.sort.bam | less -S
```

## Alignment QC

As a first QC, we can run samtools stats:

```
samtools stats alignment/SL2_CoV.sort.bam > stats/SL2_samtools_stats.txt
```

-	How many reads were mapped?
-	Which was the average length of the reads? And the maximum read length?

Additionally, NanoStat can also be used:

```
NanoStat --bam alignment/SL2_CoV.sort.bam > stats/SL2_nanostat_stats.txt
```

Now we will get the coverage per base using samtools depth.

```
samtools depth alignment/SL2_CoV.sort.bam > stats/coverage.txt
```

And look at the coverage of the SL2 sample across the entire reference genome. For that, you can start R from the command-line:

```
R
```

and then, type the following:
```
library(ggplot2)
coverage <-  read.table("stats/coverage.txt", header=FALSE, col.names = c("chrom", "pos", "cov"))
p1 <- ggplot(coverage, aes(x = pos, y = cov)) + 
      geom_line() + 
      scale_y_log10()
p1
```

Additionally, you can also look at the coverage distribution in R:

```
cov_percent <- data.frame(  "cov" = seq(1,max(coverage$cov)) 
                          , "percent" = sapply(seq(1,max(coverage$cov)), function(x) nrow(coverage[coverage$cov >= x,])/nrow(coverage)))
p2 <- ggplot(cov_percent, aes(x = cov, y = percent)) + 
      geom_line() + 
      scale_x_continuous(breaks=seq(0,max(coverage$cov), 1000)) + 
      xlab("Coverage") + 
      ylab("Percentage of bases")
p2
```

You can also add a vertical line to the previous plot intercepting the median coverage:

```
p2 + geom_vline(xintercept=median(coverage$cov), colour = "red") +
     geom_text(aes(y = .5, x = 700, label = median(coverage$cov), colour = "red"))
```

## Alignment visualisation

To inspect the alignment, we will use [Integrative Genomics Viewer](https://software.broadinstitute.org/software/igv/).

Open IGV and load the HCov-229E.fasta reference genome by selecting Genomes>Load from File or Genomes>Load from URL. The new genome will be added to the drop-down menu, and also loaded and displayed in the IGV window.

Then, load your BAM file located in alignment/SL2_CoV.sort.bam


