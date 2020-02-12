# 5. Variant calling

In this section we will cover:

* [Working directory](#workingdirectory)
* [Single Nucleotide Variant Calling](#snvcalling)
* [Structural Variant Calling](#svcalling)
* [Structural Variant Annotation](#svannotation)

You will learn to:

- Call SNVs from nanopore data and compare them to short-read sequencing SNVs
- Call SVs and interpret the results
- Make nice reports

## Working directory {#workingdirectory}

Before we start, go to your wd:

```
cd ~/Course_Materials/wd/day2
```

And make the following directories:

```
mkdir variant_calling
mkdir variant_calling/snvs
mkdir variant_calling/svs
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

*Note:* This step will take ~18min!! Feel free to grab a coffee, ask questions or continue to the next section (SV calling) in a different Terminal window - you can always come back later :-)

...

Wellcome back! When the SNV calling is done, compress and index the VCF file:

```
bgzip variant_calling/snvs/LRS_SNVs.vcf
tabix -p vcf variant_calling/snvs/LRS_SNVs.vcf.gz
```

Now you can compare the instersection between both LRS- and SRS-VCF files using bcftools:

```
bcftools isec -p isec variant_calling/snvs/LRS_snvs.vcf.gz $SRS_snvs
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

*Hint:* to count the number of variants (= numer of rows in file excluding header) you can use the command `bcftools view -H isec/0000.vcf | wc -l`

Now, explore the variants in the **IGV**. Open IGV and load the Homo_sapiens.GRCh38.dna.fasta (`$ref`) reference genome by selecting Genomes>Load from File or Genomes>Load from URL. The new genome will be added to the drop-down menu, and also loaded and displayed in the IGV window.

Then, load your BAM files `$LRS_bam` and `$SRS_bam`.

Go to positions:
- 


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

**Note: The commands above are only for the VCF run with default parameters. Do the same for the `s20` VCF.**

Now, we will compare how many SVs called in each VCF:

```
bcftools view -H variant_calling/svs/LRS_SVs.sort.vcf.gz | wc -l
```

To investigate the differences, we are going to intersect both VCF files:

```
bcftools isec -p variant_calling/svs/isec variant_calling/svs/LRS_SVs.sort.vcf.gz variant_calling/svs/LRS_SVs_s20.sort.vcf.gz
```

Do the same for the `s20` file.

- How many variants were called in both VCF files?
- Why are there variants in the `s20` VCF file that are not in the default one?
- Open both VCF files in IGV and inspect the following regions (all of them absent in `s20`)
     - 7:89890594-89901123
     - 7:90462621-90476627
- Do you think that `s20` is too strict or lenient?


## Structural Variant Annotation {#svannotation}

In order to perform annotation of the SVs from multiple sources, we will use [**AnnotSV**](https://lbgi.fr/AnnotSV).

AnnotSV annotates SVs with information about the genes (OMIM, ClinGen), regulatory elements (enhancders, promoters), pathogenicity (known from dbVar), frequencies (gnomAD, internal) and breakpoints (GC content, repeats) they overlap.

To run a basic SV annotation, we will exectue the following command:

```
AnnotSV -SVinputFile variant_calling/svs/LRS_SVs.sort.vcf.gz -SVinputInfo 1 \
     -genomeBuild GRCh38 \
     -outputDir annotation
     -outputFile LRD_SVs.anno.vcf
```

Now... Inspect the calls in *SGCE* gene:

```
grep SGCE annotation/LRD_SVs.anno.vcf.tsv
```

And visualise them in IGV.

Mutations in *SGCE* gene [[MIM:159900]](https://www.omim.org/entry/159900) have previously been associated with Myoclonic Dystonia. As you can see in the IGV plot, there are multiple deletions in this individual overlaping *SGCE*. This data comes from a patient with Myoclonic Dystonia. The deletions you just found have already been seen to be pathogenic. Congratulations, you just called variants from long-read sequencing data and identified a pathogenic one associated with disease! And actually... a very complex one! Wait for the wrap-up :)
