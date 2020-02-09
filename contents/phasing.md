#  Oxford Nanopore - day 2 - Phasing

#### [Phasing - WhatsHap](#Phasing-WA) <br>
#### [Phasing - Visualization](#Phasing-Visualization) <br>
#### [Phasing - Pedigrees](#Phasing-Pedigrees)

***

N.B.: Note these scripts use `path/to/course/` as working directory

## Phasing - WhatsHap {#Phasing-WA}

For this practical we are going to use the [WhatsHap tool](https://www.biorxiv.org/content/10.1101/085050v2.full.pdf) based on the homonym algorithm. 

[WhatsHap](https://whatshap.readthedocs.io/en/latest/index.html) is a python languaged based phasing tool. It can phase SNVs and indels, but it does not work on structural variants. WhatsHap is very easy to use: it takes as input a BAM file and a VCF file and returns a second VCF file improved with phasing information. 

BAM and VCF files don't have to be derived from the same set of reads. This is convenient bacause high-quality short reads (such us Illumina reads) can be used for variant calling and then these can be paired to long reads (such us Oxford Nanopore reads) for the phasing. Indeed, the workflow that has just been described is the reccomended one. 

However, the use of 2 different set of reads may creates some problem to the WhatHap tool, which is trying to match the sample name (`@RG` header line) of the VCF file to the BAM file, but we will see a couple of way to solve this problem.

The main WhatsHap subcommand is `phase`, which is the backbone of the entire tool. In this example we are going to phase the VCF we use as test set.

```{}
# The bamtools needs to be indexed bedore, so we do:
samtools index phasing/res/sample1_long_reads.bam

whatshap phase --indels --ignore-read-groups --reference phasing/res/reference_chr20.fa -o phasing/sample1_phased.vcf phasing/res/sample1_short_reads.vcf.gz phasing/res/sample1_long_reads.bam
```

`--ignore-read-groups` is to ask `whatshap` do overlook teh sample name information. But, because `whatshap` is not checking anymore that the samples name are teh same we have to be sure that the iput file are referring to tha same sample and that the BAM file and VCF file are referring only to one sample.

You can use the following commands to cehck the sample name in the 2 files:
```{}
samtools view -H phasing/res/sample1_long_reads.bam | grep "@RG"
# sample name isthe one starting with "SM:"

bcftools query -l phasing/res/sample1_short_reads.vcf.gz 
```

WhatsHap was able to discriminate the allele for 4601 variant as it printed out in the report (see `Found 19244 reads covering 4601 variants`). If we now look at SNVs and how they are distributed on the alleles we should that the command has divided them on allele 1 and allele 2 (for the one where this was possible):

```{}
# Before the phasing
bcftools view -H phasing/res/sample1_short_reads.vcf.gz | grep -c "0/1"
bcftools view -H phasing/res/sample1_short_reads.vcf.gz | grep -c "0|1"
bcftools view -H phasing/res/sample1_short_reads.vcf.gz | grep -c "1/0"
bcftools view -H phasing/res/sample1_short_reads.vcf.gz | grep -c "1|0"

# After the phasing
bcftools view -H phasing/sample1_phased.vcf | grep -c "0/1"
bcftools view -H phasing/sample1_phased.vcf | grep -c "0|1"
bcftools view -H phasing/sample1_phased.vcf | grep -c "1/0"
bcftools view -H phasing/sample1_phased.vcf | grep -c "1|0"
```

The different genotype separator (i.e. `|` or `/`) is because different variant caller use different separators, but both are equally good in the vcf notation. Please note that now the SNVs and indel are roughly 50% on allele 1 (a.k.a. `1|0`) and 50% on allele 2 (a.k.a `0|1`), while before the pahsing the were mostly on allele 2.


## Phasing - Visualization {#Phasing-Visualization}

### Visualise Haplotypes

It is always good practive to visually inspect our data, for this purpose we could use IGV. WhatsHap has a script to convert the VCF annotation to a GTF format that can be loaded on IGV.

```{}
whatshap stats --gtf=phasing/sample1_phased.gtf phasing/sample1_phased.vcf
```

This command will print out some stats about the haplotype blocks called in the previous command and it will also create a gtf file that we can load on IGV. Open the IGV browser and load the `phasing/sample1_phased.gtf` file you just created

Go on chr20:40000000-45000000 to look at the subset we are using.

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/HaploBlocks.png" alt="img_1" class="inline"/>

### Visualise Grouped Reads

We could also visualise haplotypes blocks and alleles directly on the reads. For this, we will need the VCF file we created and the bam file used to create it. the command is adding 2 tags to the bam file `PS` and `HP` which stasds for phase and haplotype respectively. On IGV we can use this tag to colour and order our reads.

```{}
bcftools view -O z -o phasing/sample1_phased.vcf.gz phasing/sample1_phased.vcf
bcftools index phasing/sample1_phased.vcf.gz
mv phasing/sample1_phased.vcf.gz.csi phasing/sample1_phased.vcf.gz.tbi

whatshap haplotag  --ignore-read-groups -o phasing/sample1_haplotagged.bam --reference phasing/res/reference_chr20.fa phasing/sample1_phased.vcf.gz phasing/res/sample1_long_reads.bam 

samtools index phasing/sample1_haplotagged.bam
```

Now, open the `phasing/sample1_haplotagged.bam` file on IGV.
`Color alignment by` > `tag` and insert `PS` and `Group alignment by` > `tag` and insert `PS` 

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/ColorAlignmentTag.png" alt="img_1" class="inline"/>
<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/ColorGroupTag.png" alt="img_1" class="inline"/>
<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/GroupAlignmentTag.png" alt="img_1" class="inline"/>

## Phasing - Pedigrees {#Phasing-Pedigrees}

If we have sequence information coming from our sample relative, we could use them to 
```{}
family_ID, proband_id, paternal_id, maternal_id, sex, phenotype
FAM_01	EXAMPLE_1	EXAMPLE_2	EXAMPLE_3	1	unknown
```

we need to merge the VCF file to one and also the bam file 

```{}
samtools index phasing/res/sample2_short_reads.bam
samtools index phasing/res/sample3_short_reads.bam
samtools merge -r -O BAM phasing/merged_pedigree_samples.bam phasing/res/sample1_long_reads.bam phasing/res/sample2_short_reads.bam phasing/res/sample3_short_reads.bam
samtools index phasing/merged_pedigree_samples.bam

bcftools index phasing/res/sample1_short_reads.vcf.gz
bcftools index phasing/res/sample2_short_reads.vcf.gz
bcftools index phasing/res/sample3_short_reads.vcf.gz
bcftools merge -O z -o phasing/merged_pedigree_samples.vcf.gz --merge all phasing/res/sample1_short_reads.vcf.gz phasing/res/sample2_short_reads.vcf.gz phasing/res/sample3_short_reads.vcf.gz  
```

```{}
whatshap phase --ped phasing/res/PED_ONT.txt -o phasing/pedigree_phased.vcf --reference phasing/res/reference_chr20.fa phasing/merged_pedigree_samples.vcf.gz phasing/merged_pedigree_samples.bam
```

Once we have phased the pedigree, the output will looks something like this:
```{}
# NB this VCF file has been simplified 
CHROM   POS             ID                       REF     ALT      FORMAT    SAMPLE
chr20   44705582        rs373596784              AAT     A        GT        0|1 # Alternative allele comes from the mother
chr20   40003760        rs6071894                G       T        GT        1|0 # Alternative allele comes from the father
chr20   44928864      	rs57778806;rs796352932	 CAA	   CA,C	 	  GT        1|2 # Both alleles are different from the reference genome, where 1 comes from paternal and 2 from maternal allele
```

The `paternal_allele | maternal_allele` used in whatshap is a convention started with the 1000 genome project and kept afterwards. 

Open the `phasing/pedigree_phased.vcf` file on IGV and have a look at the vcf file now. It reports the anno

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/PedigreePhase.png" alt="img_1" class="inline"/>
## Phasing - Excercise

If you made it this far and you still spare some time, there is a quick exercise you may want to do. For sample2 and sample3 we have VCF and BAM file. However, the bam files are derived from Illumina sequencing. You can try to phase this samples ans cehck the haplotype you get to the one we got for sample1

