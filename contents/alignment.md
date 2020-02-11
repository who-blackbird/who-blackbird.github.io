# 4. Alignment

In this section we will cover:

* [Background](#Background)
* [Data](#Data)
* [WorkingDirectory](#WorkingDirectory)
* [ReadsQC](#ReadsQC)
* [Alignment](#Alignment)
* [AlignmentQC](#AlignmentQC)
* [Visualisation](#Visualisation)

You will learn to:
- Take the basecalled sequences in FASTQ format and align them to a genome of reference
- Perform quality control and visualise the alignment

## Background

Coronaviruses are enveloped positive-sense (+) single-stranded RNA viruses that infect a variety of mammalian and avian hosts. Among the several coronaviruses that are pathogenic to humans, most are associated with mild clinical symptoms, with the exception of severe acute respiratory syndrome (SARS) coronavirus (SARS-CoV), a Middle East respiratory syndrome (MERS) coronavirus (MERS-Cov), and the recently identified human-infecting coronavirus, provisionally named 2019 novel coronavirus (2019-nCoV).

In this practical we will be using data from Human Coronavirus 229E (HCoV-229E), one of the first coronavirus strains being described, with a genome size of ~27,300nt. In HCov-229E-infected cells, a total of seven major viral RNAs are produced. In it's 5'-terminal region, the genome RNA contains two large ORFs, 1a and 1b, that encode the viral replicase polyproteins 1a and 1b. mRNAs 2, 4, 5, 6 and 7 are used to produce the S protein, accessory protein 4, E protein, M protein and N protein respectively.

The scheme of genomic and subgenomic RNAs produced in HCov-229E-infected cells looks like this:

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/HCov-RNAs_CoV.png" alt="rnas" class="inline"/>

## Data

The data we will be using is from Human Coronavirus 229E (HCoV-229E), one of the first coronavirus strains being described.

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

## WorkingDirectory

Now, we will set up the working directory where we will do the analysis. Open your terminal, go to 

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

## ReadsQC

There are many approaches to assess the quality of the reads. Here we will use [NanoStat](https://github.com/wdecoster/nanostat). This calculates various statistics from a long read sequencing dataset in FASTQ, BAM or albacore sequencing summary format.

For FASTQ files, it provides information for the number of reads, the read lenght and quality distrubutions (including mean and median) and the read length N50.

From your wd, run:

```
NanoStat --fastq ~/Course_Materials/nanopore_practical/data/fastq/WT_CoV.fastq.gz > stats/WT_fastq_nanostats.txt
NanoStat --fastq ~/Course_Materials/nanopore_practical/data/fastq/SL2_CoV.fastq.gz > stats/SL2_fastq_nanostats.txt
```

- What are the mean and median read length and quality for the WT and SL2 samples?

However, you may want to run your customised scripts, to answer your own questions... Here you have some ideas!

For example, to calculate how many reads we have in the FASTQ files you can also use `awk`:

```
zcat ~/Course_Materials/nanopore_practical/data/fastq/WT_CoV.fastq.gz | awk '{s++}END{print s/4}' 
zcat ~/Course_Materials/nanopore_practical/data/fastq/SL2_CoV.fastq.gz | awk '{s++}END{print s/4}'
```

As well as to obtain the read length for each read:

```
zcat ~/Course_Materials/nanopore_practical/data/fastq/WT_CoV.fastq.gz  | awk '{if(NR%4==2) print length($1)}' > stats/WT_read_length.txt
zcat ~/Course_Materials/nanopore_practical/data/fastq/SL2_CoV.fastq.gz | awk '{if(NR%4==2) print length($1)}' > stats/SL2_read_length.txt
```

Then you can look at the read length distribution. For that, you can start R from the command-line:

```
R
```

and type the following:

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

#Reorder the groups
readLength$sample <- factor(readLength$sample, levels = c('WT', 'SL2'))

#Make the plot
p1 <- ggplot(data=readLength, aes(length, fill = sample)) + 
    geom_histogram(position = 'identity', alpha = 0.6, binwidth = 10) +
    scale_x_continuous(breaks = seq(0, max(readLength$length), 1000))
p1
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

Here we will use **minimap2** to map the reads to the genome of reference. We will start with the WT sample. Then we will convert the SAM output to BAM format and sort it by mapping coordinate.

```
##Align to the ref using minimap
minimap2 -x map-ont --MD -u n -k 14 -a ~/Course_Materials/nanopore_practical/data/reference_genome/HCov-229E.fasta ~/Course_Materials/nanopore_practical/data/fastq/WT_CoV.fastq.gz > alignment/WT_CoV.sam

##Sort and output as BAM
samtools sort alignment/WT_CoV.sam -O BAM -o alignment/WT_CoV.sort.bam
```

***Alternatively***, you can run these two steps using only one command line. Below is the alignment for the SL2 sample:

```
minimap2 -x map-ont --MD -u n -k 14 -a ~/Course_Materials/nanopore_practical/data/reference_genome/HCov-229E.fasta ~/Course_Materials/nanopore_practical/data/fastq/SL2_CoV.fastq.gz | samtools sort - -O BAM -o alignment/SL2_CoV.sort.bam
```

Finally we will index the BAM files to run samtools subtools later.

```
samtools index alignment/SL2_CoV.sort.bam
samtools index alignment/WT_CoV.sort.bam
```

To visualise a BAM file:

```
samtools view alignment/SL2_CoV.sort.bam | less -S
```

## AlignmentQC

As a first QC, we can run NanoStat on the BAM files:

```
NanoStat --bam alignment/WT_CoV.sort.bam > stats/WT_bam_nanostat.txt
NanoStat --bam alignment/SL2_CoV.sort.bam > stats/SL2_bam_nanostat.txt
```
Compare the
- Mean and median aligned read length
- Number of reads aligned
- Read length N50
- Quality cutoffs
- Top highest quality and longest reads


Another option is to run samtools stats:

```
samtools stats alignment/WT_CoV.sort.bam > stats/WT_samtools_stats.txt
samtools stats alignment/SL2_CoV.sort.bam > stats/SL2_samtools_stats.txt
```

-	How many reads were mapped and unmapped? Why do you think is that?

To obtain the coverage per base, we can use `samtools depth`:

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
    scale_x_continuous(breaks = seq(0, max(coverage$pos), 1000))
p2
```

Additionally, you can also look at the coverage distribution in R:

```
percent_cov <- coverage %>% 
    group_by(sample) %>% 
    mutate(percent = (map_int(cov, ~ sum(cov >= .x)))/27317) %>% 
    select(sample, cov, percent) %>% 
    unique
    
#Make the plot
p3 <- ggplot(data = percent_cov, aes(x = cov, y = percent*100, colour = sample)) + 
    geom_line() + 
    scale_x_continuous(breaks=seq(0,max(coverage$cov), 1000)) + 
    xlab("Coverage") + 
    ylab("Percentage of bases")
p3
```

You can also add a vertical line to the previous plot intercepting with a minimum coverage of (let's say) 100x:

```
p3 + geom_vline(xintercept = 100, colour = "red")
```

## Visualisation

To inspect the alignment, we will use [Integrative Genomics Viewer](https://software.broadinstitute.org/software/igv/).

Open IGV and load the HCov-229E.fasta reference genome by selecting Genomes>Load from File or Genomes>Load from URL. The new genome will be added to the drop-down menu, and also loaded and displayed in the IGV window.

Then, load your BAM file located in alignment/SL2_CoV.sort.bam
