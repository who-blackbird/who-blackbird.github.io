# 5. Variant calling

In this section we will cover:

- [Working directory](#workingdirectory)
- [Single Nucleotide Variant Calling](#snvcalling)
- [Structural Variant Calling](#svcalling)
- [Make a Report](#report)
- [Structural Variant Annotation](#svannotation)
- [Represent cxSVs using Circos](#circos)

You will learn to:

- Call SNVs from nanopore data and compare them to short-read sequencing SNVs
- Call SVs and interpret the results
- Make a report with the results form the variant calling
- Annotate SVs with gene information
- Make a super cool Circos plot

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/vcalling_schema.png" alt="vcalling" class="inline"/>

## Working directory {#workingdirectory}

Before we start, set up your wd and go into it:

```
mkdir -p ~/Course_Materials/wd/day2
cd ~/Course_Materials/wd/day2
```

And make the following directories:

```
mkdir variant_calling
mkdir variant_calling/snvs
mkdir variant_calling/svs
mkdir annotation
mkdir circos
```

You can also define the following variables that we will use later for convienience:

```
LRS_bam=~/Course_Materials/data/variant_calling/bams/LRS_alignment.bam
SRS_bam=~/Course_Materials/data/variant_calling/bams/SRS_alignment.bam
SRS_snvs=~/Course_Materials/data/variant_calling/snvs/SRS_SNVs.vcf.gz
ref=~/Course_Materials/data/variant_calling/reference_genome/Homo_sapiens.GRCh38.dna.fasta
```

## Single Nucleotide Variant Calling {#snvcalling}

In this section we will identify SNVs from long-read sequencing data using [**Longshot**](https://github.com/pjedge/longshot).

Variants are called and stored in [VCF](http://samtools.github.io/hts-specs/VCFv4.2.pdf) format. This contains a header, and then data lines each containing information about a position in the genome.

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/vcf.png" alt="img_3" class="inline"/>

Call SNVs in the nanopore aligment example as below:

```
longshot --bam $LRS_bam --ref $ref --out variant_calling/snvs/LRS_SNVs.vcf
```

When the SNV calling is done, compress and index the VCF file:

```
bgzip variant_calling/snvs/LRS_SNVs.vcf
tabix -p vcf variant_calling/snvs/LRS_SNVs.vcf.gz
```

Now you can compare the intersection between both LRS- and SRS-VCF files using bcftools:

```
bcftools isec -p variant_calling/snvs/isec variant_calling/snvs/LRS_SNVs.vcf.gz $SRS_snvs
```

This will create a folder named isec with the following files:

```
isec/0000.vcf   for records private to the long_reads_VCF
isec/0001.vcf   for records private to the short_reads_VCF
isec/0002.vcf   for records from long_reads_VCF shared by both
isec/0003.vcf   for records from short_reads_VCF shared by both
```

- How many SNVs have been called by both technologies?
- How many SNVs have been missed by short and/or long read sequencing?

_Hint:_ to count the number of variants (= number of rows in file excluding header) you can use the next command:

```
for f in variant_calling/snvs/isec/000*.vcf
do
    echo $f `bcftools view -H $f | wc -l`
done
```

Note that most of the variant that are private to short-reads are indels, since longshot does not call them.

To get the number of SNPs called by short-reads, and missed by long-reads, run:

```
bcftools view -HV indels variant_calling/snvs/isec/0001.vcf | wc -l
```

And to get how many of those are only PASS variants:

```
bcftools view -f PASS -HV indels variant_calling/snvs/isec/0001.vcf | wc -l
```

Let's explore them in the **IGV**! Open IGV and load the Homo_sapiens.GRCh38.dna.fasta (`$ref`) reference genome by selecting Genomes>Load from File or Genomes>Load from URL. The new genome will be added to the drop-down menu, and also loaded and displayed in the IGV window.

Then, load your BAM files `$LRS_bam` and `$SRS_bam`.

Feel free to look at any of the variants. We also suggest:
- 7:88782662
- 7:88867497
- 7:88874102
- 7:89194947


## Structural Variant Calling {#svcalling}

Algorithms for calling SVs from long-read sequencing data include:

- [Sniffles](http://github.com/fritzsedlazeck/Sniffles): best used with minimap2 or NGMLR.
- [NanoSV](http://github.com/philres/ngmlr): best used with LAST.

Here, we will use **Sniffles** for calling structural variants.

```
sniffles -m $LRS_bam -v variant_calling/svs/LRS_SVs.vcf
```

If you want to look at high quality SVs, you can change the -s parameter to 20, where s is the minimum number of reads that support a SV (by default is 10).

```
sniffles -m $LRS_bam -v variant_calling/svs/LRS_SVs_s20.vcf -s 20
```

The information that is provided in sniffles’s output can be found in:
[http://github.com/fritzsedlazeck/Sniffles/wiki/Output](http://github.com/fritzsedlazeck/Sniffles/wiki/Output)

Sort the VCF files (for that you'd need to compress them first) and index them:

```
bgzip variant_calling/svs/LRS_SVs.vcf
vcf-sort variant_calling/svs/LRS_SVs.vcf.gz | bgzip -c > variant_calling/svs/LRS_SVs.sort.vcf.gz
tabix -p vcf variant_calling/svs/LRS_SVs.sort.vcf.gz
```

Do the same for the `s20` VCF

```
bgzip variant_calling/svs/LRS_SVs_s20.vcf
vcf-sort variant_calling/svs/LRS_SVs_s20.vcf.gz | bgzip -c > variant_calling/svs/LRS_SVs_s20.sort.vcf.gz
tabix -p vcf variant_calling/svs/LRS_SVs_s20.sort.vcf.gz
```

Now, we will compare how many SVs were called in each VCF:

```
bcftools view -H variant_calling/svs/LRS_SVs.sort.vcf.gz | wc -l
bcftools view -H variant_calling/svs/LRS_SVs_s20.sort.vcf.gz | wc -l
```

To investigate the differences, we are going to intersect both VCF files:

```
bcftools isec -p variant_calling/svs/isec variant_calling/svs/LRS_SVs.sort.vcf.gz variant_calling/svs/LRS_SVs_s20.sort.vcf.gz
```

Run the next command to know how many variants have been identified in each VCF file:

```
for f in variant_calling/svs/isec/000*.vcf
do
    echo $f `bcftools view -H $f | wc -l`
done
```

- How many variants were called in both VCF files?
- Why are there variants in the `s20` VCF file that are not in the default one?
- Open both VCF files in IGV and inspect the following regions (all of them absent in `s20`)
  - 7:89890594-89901123
  - 7:89969558-89969874
  - 10:52936870-53686004
- Do you think that `s20` is too strict or lenient?

## Make a Report {#report}

We will now make a report, like the one you made in the Quality Control section, but for the SVs you just called.

For that, we will use an already existing [tutorial](https://github.com/nanoporetech/ont_tutorial_sv).

Open the `config.yaml` file:

```
cat ~/Course_Materials/data/variant_calling/ont_tutorial_sv/config.yaml
```

That should look like this:

```
---
# this config.yaml is passed to Snakefile in pipeline-structural-variation subfolder.
# Snakemake is run from this pipeline-structural-variation folder; it is necessary to
# pass an appropriate path to the input-files (the ../ prefix is sufficient for this demo)

# FASTA file containing the reference genome
reference_fasta: "~/Course_Materials/data/variant_calling/reference_genome/Homo_sapiens.GRCh38.dna.fasta"
# Sample name
sample_name: "SAMPLE"

##################################
## Tutorial specific parameters ##
##################################
tutorialText: FALSE
biocGenome: "hg38"
GeneIdMappings: "org.Hs.eg.db"
GenomeAnnotation: "TxDb.Hsapiens.UCSC.hg38.knownGene"
RepeatDB: "~/Course_Materials/data/variant_calling/repeatmasker/hg38.fa.out.gz"
bamFile: "~/Course_Materials/data/variant_calling/bams/LRS_alignment.bam"
vcfFile: "~/Course_Materials/wd/day2/variant_calling/svs/LRS_SVs.vcf.gz"
```

This is the config file that points to the path of the files that it requires.

To make the report, we are interested in the `ont_tutorial_sv.Rmd` file. 

Open Rstudio and open it.

Press the “Knitr” button at the top of the page and wait for magic!

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/qc.rhelp.png" alt="help" class="inline"/>

Again... if magic does not happen, R should have printed a .html file - open it by pasting the below in the terminal!

```
open ~/Course_Materials/data/variant_calling/ont_tutorial_sv/ont_tutorial_sv.html
```

## Structural Variant Annotation {#svannotation}

In order to perform annotation of the SVs from multiple sources, we will use [**AnnotSV**](https://lbgi.fr/AnnotSV).

AnnotSV annotates SVs with information about the genes (OMIM, ClinGen), regulatory elements (enhancers, promoters), pathogenicity (known from dbVar), frequencies (gnomAD, internal) and breakpoints (GC content, repeats) they overlap.

To run a basic SV annotation, we will execute the following command:

```
AnnotSV -SVinputFile variant_calling/svs/LRS_SVs.sort.vcf.gz -SVinputInfo 1 \
     -genomeBuild GRCh38 \
     -outputDir annotation \
     -typeOfAnnotation full
```

Now... Inspect the calls in _PCDH15_ gene:

```
grep PCDH15 annotation/LRD_SVs.anno.vcf.tsv
```

And visualise them in IGV.

Next, you can also look at the calls overlapping the _PCDH15_ gene.

```
grep PCDH15 annotation/LRD_SVs.anno.vcf.tsv
```

What do you think is happening there? 

...


Yes! This is a complex SV!


## Represent cxSVs using Circos {#circos}

The previous cxSV involves multiple breakpoints across different chromosomes. To characterise it, you are going to represent it using [**Circos**](http://circos.ca).

Circos needs the following input files:
- Scaffold
- Links
- Config file
- Coverage

To save you some time, we already prepared these files and put them under ```~/Course_Materials/data/circos```.

Open the files to see the information they contain.

Now, you can make the plot using circos software:

```
cd circos
circos -conf ~/Course_Materials/data/circos/circos.conf
```

Yes, it is true - making circos plots is not as difficult as it looks like!! Actually, it is quite cool :-)

You can now visualise the plot you just made:

```
xdg-open circos.png
```

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/circos.png" alt="circos" class="inline"/>
