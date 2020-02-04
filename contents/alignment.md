# 4. Alignment

## Working directory

Open your terminal, and go to your working directory, 

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

## Variables

First, define your variables:

```
fastq=~/Course_Materials/nanopore_practical/data/Viehweger_Jena/minion/fastq/2017-09-29_coronavirus_SL2.rna.fastq.gz
ref=~/Course_Materials/nanopore_practical/data/reference_genome/HCov-229E.fasta
```

## Data

The data we will be using is from [NA12878](http://github.com/nanopore-wgs-consortium/NA12878/blob/master/Genome.md) human genome reference standard on the Oxford Nanopore MinION using 1D ligation kits (450 bp/s) and R9.4 chemistry (FLO-MIN106).

We have already prepared a subset of specific regions of NA12878 genome in a FASTQ file. [FASTQ](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2847217) format is a text-based format for storing both a biological sequence and its corresponding quality scores.

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/fastq.png" alt="img_1" class="inline"/>

A FASTQ file normally uses four lines per sequence: 
 1) Begins with a ‘@’ and is followed by a sequence identifier 
 2) Is the raw sequence letters
 3) Begins with a ‘+’ character 
 4) Encodes the quality values for the sequence in Line 2

You can visualize the FASTQ file typing:

```
less -S $fastq
```

## Reads QC

First we will calculate how many reads we have:

```
awk '{s++}END{print s/4}' $fastq
```

You can then use awk to obtain the read length for each read:

```
awk '{if(NR%4==2) print length($1)}' $fastq > stats/read_length.txt
```

And look at the read length distribution. For that, you can start R from the command-line:

```
R
```

and then, type the following:

```
library(ggplot2)
readLength <- read.table("stats/read_length.txt", header=FALSE, col.names = "length")
head(readLength)
ggplot(data=readLength, aes(length)) + geom_histogram()
```

For quitting R, just type:

```
quit()
```

## Alignment

The standard format for aligned sequence data is [SAM](http://samtools.github.io/hts-specs/SAMv1.pdf) (Sequence Alignment Map). 

SAM files have a header that contains information on alignment and contigs used, and the aligned reads:

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/sam.jpg" alt="img_2" class="inline"/>

But because SAM files can be large, they are usually stored in the compressed version of them, [BAM](http://samtools.github.io/hts-specs/SAMv1.pdf) files.

Multiple algorithms have been developed to align long reads to a genome of reference. Some examples are:
-	Graphmap: [http://github.com/isovic/graphmap](http://github.com/isovic/graphmap)
-	bwa mem -x l ont2d: [http://github.com/lh3/bwa](http://github.com/lh3/bwa)
-	LAST: [http://last.cbrc.jp](http://last.cbrc.jp)
-	NGMLR: [http://github.com/philres/ngmlr](http://github.com/philres/ngmlr)
-	minimap2: [http://github.com/lh3/minimap2](http://github.com/lh3/minimap2)

Here we will use minimap2 to map the reads to the genome of reference, and convert the SAM output to BAM format.

```
minimap2 -x map-ont --MD -a $ref $fastq > alignment/SL2_CoV.sam
samtools view alignment/SL2_CoV.sam -O BAM -o alignment/SL2_CoV.bam
```

Then, we will sort it by mapping coordinate and save it as BAM.

```
samtools sort alignment/SL2_CoV.bam > alignment/SL2_CoV.sort.bam
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


