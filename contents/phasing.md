#  Oxford Nanopore - day 2 - Phasing

#### [Phasing - WhatsHap](#Phasing-WA) 
#### [Phasing - Visualization](#Phasing-Visualization) 
#### [Phasing - Pedigrees](#Phasing-Pedigrees)

***

N.B.: Note these scripts use `path/to/course/` as working directory

## Phasing - WhatsHap {#Phasing-WA}

For this practical we are going to use the [WhatsHap tool](https://www.biorxiv.org/content/10.1101/085050v2.full.pdf), which uses python and is based on the homonym algorithm. To download [WhatsHap](https://whatshap.readthedocs.io/en/latest/index.html). 

WhatsHap can phase SNVs and indels, but it does not work on structural variants. It is very easy to use: it takes as input a BAM file and a VCF file and returns a second VCF file improved with the phasing information. BAM and VCF files don't have to be derived from the same set of reads. This is convenient bacause one can use high-quality short reads (such us Illumina reads) to call the variants and then use long reads (such us Oxford Nanopore reads) to phase them. Indeed, this workflow that has just been described is the raccomended one. Please note that using 2 different sets of reads may create some problems to the WhatsHap algortihm, which is trying to match the sample name (`@RG` header line) of the BAM file to the VCF file. It tries to match them because this tool can process VCF and BAM containing different samples in the same file. Nonetheless, there are a couple of ways to work around this problem.

The main WhatsHap subcommand is `phase`, which is the backbone of the entire tool, it uses reads containing 2 or more variants to assign them to one allele or the other. In this example we are going to phase a test VCF.

```{}
mkdir -p Course_Materials/wd/day2/phasing/
IN="Course_Materials/data/phasing"
OUT="Course_Materials/wd/day2/phasing"
```


```{}
# The bamtools needs to be indexed bedore, so we do:
samtools index ${IN}/res/sample1_long_reads.bam

whatshap phase --indels --ignore-read-groups --reference ${IN}/res/reference_chr20.fa -o ${OUT}/sample1_phased.vcf ${IN}/res/sample1_short_reads.vcf.gz ${IN}/res/sample1_long_reads.bam
```

`--ignore-read-groups` asks `whatshap` to overlook the sample name information. But, because of this we have to be sure that the input files refer to the same sample.

You can use the following commands to check the sample names in the 2 files:
```{}
samtools view -H ${IN}/res/sample1_long_reads.bam | grep "@RG"
# sample name isthe one starting with "SM:"

bcftools query -l ${IN}/res/sample1_short_reads.vcf.gz 
```

WhatsHap was able to discriminate the allele for 4601 variants as it printed out in the report (see `Found 19244 reads covering 4601 variants`). If we now look at SNVs and how they are distributed on the alleles, we should see that the command has divided them into allele 1 and allele 2 (when this was possible):

```{}
# Before the phasing
bcftools view -H ${IN}/res/sample1_short_reads.vcf.gz | grep -c "0/1"
bcftools view -H ${IN}/res/sample1_short_reads.vcf.gz | grep -c "0|1"
bcftools view -H ${IN}/res/sample1_short_reads.vcf.gz | grep -c "1/0"
bcftools view -H ${IN}/res/sample1_short_reads.vcf.gz | grep -c "1|0"

# After the phasing
bcftools view -H ${OUT}/sample1_phased.vcf | grep -c "0/1"
bcftools view -H ${OUT}/sample1_phased.vcf | grep -c "0|1"
bcftools view -H ${OUT}/sample1_phased.vcf | grep -c "1/0"
bcftools view -H ${OUT}/sample1_phased.vcf | grep -c "1|0"
```

There are 2 different genotype separators (i.e. `|` or `/`) because different variant callers use different separators, but both are equally good in the VCF grammar. Please note that now the SNVs and indels are roughly 50% on allele 1 (a.k.a. `1|0`) and 50% on allele 2 (a.k.a `0|1`), while before the pahsing the were mostly on allele 2.

## Phasing - Visualization {#Phasing-Visualization}

### Visualise Haplotypes

It is always good practice to visually inspect our data, for this purpose we could use IGV. WhatsHap has a script to convert the VCF annotation to a GTF format that can be loaded on IGV.

```{}
whatshap stats --gtf=${OUT}/sample1_phased.gtf ${OUT}/sample1_phased.vcf
```

This command will print out some stats about the haplotype blocks called in the previous command and it will also create a GTF file that we can load on IGV. Open the IGV browser and load the `${OUT}/sample1_phased.gtf` file you just created.

Go on chr20:40000000-45000000 to look at the subset we are using.

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/HaploBlocks.png" alt="img_1" class="inline"/>

### Visualise Grouped Reads

We could also visualise haplotype blocks and alleles directly on the reads. For this, we will need the VCF file we created and the BAM file used to create it. The command adds 2 tags to the bam file `PS` and `HP` which stand for `phase` and `haplotype` respectively. On IGV we can use this tag to colour and order our reads.

```{}
bcftools view -O z -o ${OUT}/sample1_phased.vcf.gz ${OUT}/sample1_phased.vcf
bcftools index ${OUT}/sample1_phased.vcf.gz
mv ${OUT}/sample1_phased.vcf.gz.csi ${OUT}/sample1_phased.vcf.gz.tbi

whatshap haplotag  --ignore-read-groups -o ${OUT}/sample1_haplotagged.bam --reference ${IN}/res/reference_chr20.fa ${OUT}/sample1_phased.vcf.gz ${IN}/res/sample1_long_reads.bam 

samtools index ${OUT}/sample1_haplotagged.bam
```

Now, open the `${OUT}/sample1_haplotagged.bam` file on IGV.
`Color alignment by` > `tag` and insert `PS` and `Group alignment by` > `tag` and insert `PS` 

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/ColorAlignmentTag.png" alt="img_1" class="inline"/>
<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/ColorGroupTag.png" alt="img_1" class="inline"/>
<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/GroupAlignmentTag.png" alt="img_1" class="inline"/>

## Phasing - Pedigrees {#Phasing-Pedigrees}

If we have sequencing data coming from our sample's parents, we could use them to phase the VCF file according to the paternal and maternal alleles. To do this, we need BAM and VCF files of the trio and an additional file called PED. PED files have this structure:

```{}
family_ID, proband_id, paternal_id, maternal_id, sex, phenotype
```

To use this `WhatsHap` function we need to merge the VCF files into one as well as the bam files. 

```{}
# Merge the VCF files
bcftools index ${IN}/res/sample1_short_reads.vcf.gz
bcftools index ${IN}/res/sample2_short_reads.vcf.gz
bcftools index ${IN}/res/sample3_short_reads.vcf.gz
bcftools merge -O z -o ${OUT}/merged_pedigree_samples.vcf.gz --merge all ${IN}/res/sample1_short_reads.vcf.gz ${IN}/res/sample2_short_reads.vcf.gz ${IN}/res/sample3_short_reads.vcf.gz  

#Merge the BAM files
samtools index ${IN}/res/sample2_short_reads.bam
samtools index ${IN}/res/sample3_short_reads.bam
samtools merge -r -O BAM ${OUT}/merged_pedigree_samples.bam ${IN}/res/sample1_long_reads.bam ${IN}/res/sample2_short_reads.bam ${IN}/res/sample3_short_reads.bam
samtools index ${OUT}/merged_pedigree_samples.bam
```

After we created the new input files, we can run the `whatshap phase` command:

```{}
whatshap phase --ped ${IN}/res/PED_ONT.txt -o ${OUT}/pedigree_phased.vcf --reference ${IN}/res/reference_chr20.fa ${IN}/merged_pedigree_samples.vcf.gz ${IN}/merged_pedigree_samples.bam
```

Once we have phased the pedigree, the output will look something like this:

```{}
# NB this VCF file has been simplified 
CHROM   POS             ID                       REF     ALT      FORMAT    SAMPLE
chr20   40003760        rs6071894                G       T        GT        1|0 # Alternative allele comes from the father
chr20   44705582        rs373596784              AAT     A        GT        0|1 # Alternative allele comes from the mother
chr20   44928864      	rs57778806;rs796352932	 CAA	   CA,C	 	  GT        1|2 # Both alleles are different from the reference genome, where 1 comes from paternal and 2 from maternal allele
```

The `paternal_allele | maternal_allele` annotation used in WhatsHap is a convention started with the 1000 genomes project and kept afterwards. 

Open the `phasing/pedigree_phased.vcf` file on IGV and have a look at the VCF file. It reports the variants of the trio. 

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/PedigreePhase.png" alt="img_1" class="inline"/>

## Phasing - Excercise

If you made it this far and you still have time and stamina, there is a quick exercise you may want to do. For sample2 and sample3 we have VCF and BAM files. However, the BAM files are derived from Illumina sequencing. You can try to phase this samples and compare the haplotype you get to the one we got for sample1. 

