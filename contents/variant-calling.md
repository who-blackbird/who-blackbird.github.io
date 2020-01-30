# 5. Variant calling

## SNV/indel calling

Variants are called and stored in [VCF](http://samtools.github.io/hts-specs/VCFv4.2.pdf) format. This contains a header, and then data lines each containing information about a position in the genome.

<img src="//raw.githubusercontent.com/ouwehand-lab/ouwehand-lab.github.io/master/images/vcf.png" alt="img_3" class="inline"/>

## SV calling

Currently, there are different algorithms for calling SVs from long-read sequencing data, including:
-	[Sniffles](http://github.com/fritzsedlazeck/Sniffles): best used with NGMLR. 
-	[NanoSV](http://github.com/philres/ngmlr): best used with LAST.

Since we used minimap2 for the alignment, now we will use sniffles for calling structural variants.

```
sniffles -m alignment/NA12878.ROI.sort.bam -v variant_calling/NA12878.ROI.vcf
```

If you want to look at high quality SVs, you can change the -s parameter to 20, where s is the minimum number of reads that support a SV (by default is 10).

```
sniffles -m alignment/NA12878.ROI.sort.bam -v variant_calling/NA12878.ROI.s20.vcf -s 20
```

The information that is provided in sniffles’s output can be found in:
[http://github.com/fritzsedlazeck/Sniffles/wiki/Output](http://github.com/fritzsedlazeck/Sniffles/wiki/Output)

To know how many SVs have been called, we will run:

```
bcftools view -H variant_calling/NA12878.ROI.vcf | wc -l
```

Finally, you can convert the VCF to a tab format:

```
../scripts/vcf2tab.py variant_calling/NA12878.ROI.vcf
```

and inspect the deletions in IGV.

-	How many deletions are real?
-	How many SVs breakpoint junctions are within repetitive sequences?
     - For that, you would need to load Repeatmasker from server (File > Load from server > Annotations > Variation and Repeats > Repeat Masker)
