# 5. Variant calling

In this section we will cover:

* [SNV Calling](#snvcalling)
* [SV Calling](#svcalling)
* [SV Annotation](#svannotation)
* [Analysis report](#report)

You will learn to:

- Call SNVs from nanopore data and compare them to short-read sequencing SNVs
- Call SVs and interpret the results
- Make nice reports

## Workind directory

Before we start, go to your wd:

```
cd ~/Course_Materials/nanopore_practical/wd
```

And make the following directories:

```
mkdir SNVs
mkdir SVs
```

You can also define the following variables that we will use later for convienience:

```
LRS_bam=~/Course_Materials/nanopore_practical/data/day2/alignment/LRS_alignment.bam
SRS_bam=~/Course_Materials/nanopore_practical/data/day2/alignment/SRS_alignment.bam
SRS_snvs=~/Course_Materials/nanopore_practical/data/day2/SNVs/SRS_SNVs.vcf.gz
ref=~/Course_Materials/nanopore_practical/data/day2/reference_genome/Homo_sapiens.GRCh38.dna.fasta
```

## SNV Calling {#snvcalling}

In this section we will identify SNVs from long-read sequencing data using [**Longshot**](https://github.com/pjedge/longshot).

Variants are called and stored in [VCF](http://samtools.github.io/hts-specs/VCFv4.2.pdf) format. This contains a header, and then data lines each containing information about a position in the genome.

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/vcf.png" alt="img_3" class="inline"/>

Call SNVs in the nanopore aligment example as below:

```
longshot --bam $LRS_bam --ref $ref --out SNVs/LRS_SNVs.vcf
```

*Note:* This step will take ~18min!! Feel free to grab a coffee, ask questions or continue to the [next section](#svcalling) in a different Terminal window - you can always come back later :-)

...

Wellcome back! When the SNV callins is done, compress and index the VCF file:

```
bgzip SNVs/LRS_SNVs.vcf
tabix -p vcf SNVs/LRS_SNVs.vcf.gz
```

Now you can compare the instersection between both using bcftools:

```
bcftools isec -p isec SNVs/LRS_snvs.vcf.gz $SRS_snvs
```

This will creat a folder named isec with the following files:

```
isec/0000.vcf   for records private to the long_reads_VCF
isec/0001.vcf   for records private to the short_reads_VCF
isec/0002.vcf   for records from long_reads_VCF shared by both
isec/0003.vcf   for records from short_reads_VCF shared by both
```
- How many SNVs have been called by both technologies?
- How many SNVs have been missed by short and/or long read sequencing?

Now, explore the variants in the **IGV**. Open IGV and load the Homo_sapiens.GRCh38.dna.fasta (`$ref`) reference genome by selecting Genomes>Load from File or Genomes>Load from URL. The new genome will be added to the drop-down menu, and also loaded and displayed in the IGV window.

Then, load your BAM files `$LRS_bam` and `$SRS_bam`.

Go to position (PENDING)


## SV Calling {svcalling}

Algorithms for calling SVs from long-read sequencing data include:
-	[Sniffles](http://github.com/fritzsedlazeck/Sniffles): best used with minimap2 or NGMLR. 
-	[NanoSV](http://github.com/philres/ngmlr): best used with LAST.

Here, we will use sniffles for calling structural variants.

```
sniffles -m $LRS_bam -v SVs/LRS_SVs.vcf
```

If you want to look at high quality SVs, you can change the -s parameter to 20, where s is the minimum number of reads that support a SV (by default is 10).

```
sniffles -m $LRS_bam -v SVs/LRS_SVs.s20.vcf -s 20
```

The information that is provided in sniffles’s output can be found in:
[http://github.com/fritzsedlazeck/Sniffles/wiki/Output](http://github.com/fritzsedlazeck/Sniffles/wiki/Output)

Sort the VCF files (for that you'd need to compress them first) and index them:

```
bgzip SVs/LRS_SVs.vcf
vcf-sort SVs/LRS_SVs.vcf.gz | bgzip -c > SVs/LRS_SVs.sort.vcf.gz
tabix -p vcf SVs/LRS_SVs.sort.vcf.gz
```

*Note:* The commands above are only for the VCF run with default parameters. Do the same for the `s20` VCF.

Now, we will compare how many SVs called in each VCF:

```
bcftools view -H SVs/LRS_SVs.sort.vcf.gz | wc -l
```

To investigate the differences, we are going to intersect both VCF files:

```
bcftools isec -p SVs/isec SVs/LRS_SVs.sort.vcf.gz SVs/LRS_SVs.s20.sort.vcf.gz
```

Do the same for the `s20` file.

- How many variants were called in both VCF files?
- Why are there variants in the `s20` VCF file that are not in the default one?
- Open both VCF files in IGV and inspect the following regions (all of them absent in `s20`)
-- 7:89890594-89901123
-- 7:90462621-90476627
- Do you think that `s20` is too strict or lenient?


## SV Annotation {svannotation}

To annotate the VCF, 

PENDING


You can also convert the VCF to a tab format:

```
~/Course_Materials/nanopore_practical/scripts/vcf2tab.py SVs/LSR_SVs.vcf
```

and inspect the deletions in IGV.

-	How many deletions are real?
-	How many SVs breakpoint junctions are within repetitive sequences?
     - For that, you would need to load Repeatmasker from server (File > Load from server > Annotations > Variation and Repeats > Repeat Masker)

## Analysis report {#report}

The tutorial document can also be prepared from the Linux command line with the following command.

```
R --slave -e 'rmarkdown::render("ont_tutorial_sv.Rmd", "html_document")'
```

## Explore structural variation using IGV

## Make representation of complex SVs using circos

## Snakemake

## References

