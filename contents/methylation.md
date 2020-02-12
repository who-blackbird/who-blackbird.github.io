#  7. Methylation

In this section we will cover:

* [Methylation - DeepSignal](#Methylation-DeepSignal)
* [Methylation - Nanopolish](#Methylation-Nanopolish)
* [Methylation - Comparison](#Methylation-Comparison) 
* [Methylation - Visualisation](#Methylation-Visualization)

You will learn to:
- Call 5mC modification from Nanopore sequencing data
- Compare differenc methods to call 5mC
- Visualize the 5mC on IGV

***
N.B.: Note these scripts use `path/to/Course_Materials/` as working directory (wd). To check what directory you are using as `wd` do `pwd`. If the output is different from `path/to/Course_Materials` please go to `Course_Materials/` with the `cd` command.

Also here, we create the output directory and then we set the two variables: `IN` and `OUT.`
```{}
mkdir -p wd/day2/methylation/
IN="data/methylation"
OUT="wd/day2/methylation"
```

## Methylation - DeepSignal {#Methylation-DeepSignal}

Deepsignal [(doi:10.1093/bioinformatics/btz276)](https://doi.org/10.1093/bioinformatics/btz276) is a Machine learning (ML) approach to identify changes in the resistance measured by the ONT that may be identified as a 5-MethylCytosin (5mC). As all ML approaches, DeepSignal requires a model. The model has been trained by the authors and you can download it with the software. DeepSignal takes ~ 24 hours to run on 4 GPU for a human genome 30x coverage.

Deepsignal is a quite straightforward tool. However, it needs:
+ Fast5 files to work (they contain the resistance information required to call 5mC).
+ if fast5 files have a Multi-fast5 format, they need to be converted to single fast5 (you can use the command [```multi_to_single_fast5```](https://github.com/nanoporetech/ont_fast5_api)).
+ It needs the ML model. One model has been trained by the Author. Alternatively, one can train a new model, especially if one has to investigate DNA/RNA modifications that are different from 5mC. To train a new model, one can use the command ```deepsignal train```. In this case, one needs a _training set_ and a _validation set_ (N.B. They need to be 2 different data sets).

#### Call 5mC modification

```{}
deepsignal call_mods --input_path ${IN}/res/fast5_files/ --model_path ${IN}/res/model.CpG.R9.4_1D.human_hx1.bn17.sn360/bn_17.sn_360.epoch_7.ckpt --result_file ${OUT}/Met_deepsignal.tsv --reference_path ${IN}/res/reference.fasta --nproc 4 --is_gpu no
```

This command calls Cytosin modifications straight from the Fast5 data. The output should be something like:

chromosome | position | strand | pos_in_strand | readname | prob_0 | prob_1 | called_label | k_mer 

```{}
chr20	5345797	+	5345797	0bd4fa76-db91-4355-9dff-7acabda704cb	t	0.87934494	0.12065502	0	AAAAAAATCGAGAGTTG
chr20	5345957	+	5345957	0bd4fa76-db91-4355-9dff-7acabda704cb	t	0.6589491	0.3410509	0	TCTTGATCCGTTATAGC
chr20	5346125	+	5346125	0bd4fa76-db91-4355-9dff-7acabda704cb	t	0.073313385	0.92668664	1	TCATTCATCGACACACC
chr20	5346181	+	5346181	0bd4fa76-db91-4355-9dff-7acabda704cb	t	0.86940426	0.1305958	0	TGTTTCACCGTGATACA
chr20	5346431	+	5346431	0bd4fa76-db91-4355-9dff-7acabda704cb	t	0.065592036	0.93440795	1	GCTGTTTCCGCTGTCTG
chr20	5346554	+	5346554	0bd4fa76-db91-4355-9dff-7acabda704cb	t	0.8229035	0.17709653	0	CCCTCTTTCGTGCTCTG
```

The deepsignal default settings have 0 as Unmethylated and 1 as Methylated

#### Call 5mC modification frequency

```{}
python3 ${IN}/scripts/call_modification_frequency.py --input_path ${OUT}/Met_deepsignal.tsv --result_file ${OUT}/Met_frequency_deepsignal.tsv --prob_cf 0.6
```
The script to call modification frequencies comes from the [deepsignal utils](https://github.com/bioinfomaticsCSU/deepsignal/tree/master/scripts)

The frequency script output should look something like:

chromosome | position | strand | pos_in_strand | prob_0_sum | prob_1_sum | count_modified | count_unmodified | coverage | modification_frequency | k_mer

```{}
chr20	4991162	-	59453004	0.060	0.940	1	0	1	1.000	GGTGCCACCGCACTCCA
chr20	4991622	-	59452544	0.086	0.914	1	0	1	1.000	TCCTGCCTCGGCCTCCT
chr20	4991653	-	59452513	0.132	0.868	1	0	1	1.000	GCAACCTCCGCCTCCTG
chr20	4992058	-	59452108	0.094	0.906	1	0	1	1.000	CGCGATCTCGGCTCACT
chr20	4992064	-	59452102	0.056	0.944	1	0	1	1.000	TAATGGCGCGATCTCGG
```

## Methylation - Nanopolish {#Methylation-Nanopolish}

[Nanopolish](https://nanopolish.readthedocs.io/en/latest/index.html) is a software to analyse signal-level Oxford Nanopore data. It was originally designed to improve consensus sequence and polish the basecall, but today it provides a series of other tools. Nanopolish needs to access the signal-level and it uses hidden Markov model (HMM) to take into account the sorrounding resistence signal when calling the base and its modification.

To look at the raw data in an efficient way, Nanopolish needs to create an index to link the read IDs (fastq) with their signal-level data in the fast5 files. We already learned how to basecall the raw data, so we will skip this step and we go straight to the indexing:

#### Index fast5 files
We use the `index` command to link fast5 files to the fastq sequences: 
```{}
nanopolish index -d ${IN}/res/fast5_files/ ${IN}/res/alignment_output.fastq
```
We need to align the reads to the reference genome. We can use mimimap2 for this and pipe the results to samtools to sort them. Ultimately, we need to index the bam file we created:
#### Align and sort fastq files
```{}
# Align the reads the the GRCh38 subset
minimap2 -a -x map-ont ${IN}/res/reference.fasta ${IN}/res/alignment_output.fastq | samtools sort -T tmp -o ${IN}/res/alignment_output_sorted.bam

# Index the sorted output
samtools index ${IN}/res/alignment_output_sorted.bam
```
Once the preprocessing is finished, we have all the required files to start the 5mC modification call. We use the command `call-methylation` to annotate the modification.
#### Call 5mC modification
```{}
nanopolish call-methylation -t 8 -r ${IN}/res/alignment_output.fastq -b ${IN}/res/alignment_output_sorted.bam -g ${IN}/res/reference.fasta -w "chr20:5,000,000-10,000,000" > ${OUT}/Met_nanopolish.tsv
```
__NB.: This command will take about 10 minutes, so it is a good time to have a break and stretch a bit your legs.__ If the script is taking too long there is a hidden file that you can use to keep going. Do:
```{}
mv ${IN}/.Met_nanopolish.tsv ${OUT}/Met_nanopolish.tsv
```

Tip:
`-t` specifies the number of threads to use for this command. This value can be optimised on the number of cores in your CPU. To check how many threads you have, on linux, you can do:
```{}
lscpu | egrep 'Model name|Socket|Thread|NUMA|CPU\(s\)'
```
This command will print CPUs and threads. There is not an optimal number of threads, it really depends on the process you are running and it always needs a bit of optimization.


#### Call 5mC modification frequency
Now we can use the modification file to call the 5mC modification frequencies:
```{}
${IN}/scripts/nanopolish/scripts/calculate_methylation_frequency.py ${OUT}/Met_nanopolish.tsv > ${OUT}/Met_frequency_nanopolish.tsv
```

## Methylation - Comparison {#Methylation-Comparison}

#### Compare the cytosines that have been identified to be modified or not modified

Bisulphyte is nowadays commonly used to detect 5mC genome-wide. We can compare the results we got from the Nanopore methylation to the bisulphite doing. We could compare the score the softwares gave to the 5mC and the regions where they identified the modification. To compare the scores we could do:
```{}
${IN}/scripts/compare_met_score.R -d ${OUT}/Met_deepsignal.tsv -n ${OUT}/Met_frequency_nanopolish.tsv -b ${IN}/res/bisulfite.ENCFF835NTC.example.tsv -o ${OUT}/methylation_score_comparison
```
Please, feel free to open the script and look at what is it doing. This command prints out a table (i.e. `methylation_score_comparison.tsv`) and a pdf (i.e. `methylation_score_comparison.pdf`) comparing the methulation scores. The table list the scores given by every software/method. It should looks like:
```{}
Index            Deepsignal_frequency  Bisulphite_frequency  Nanopolish_frequency
chr20:5017505    0.88322896	           0.75	                 0.8
chr20:5017536	   0.46293265	           0.6	                 0.857
chr20:5017607	   0.9334746	           0.8	                 0.818
chr20:5017628  	 0.9180513	           1	                   0.818
chr20:5017662	   0.15342923	           0.78	                 0.462
chr20:5017701	   0.92652583	           0.62	                 0.778
chr20:5017906	   0.86876917	           1	                   0.846
chr20:5017998	   0.28020185	           0.29	                 0.286
```
and the pdf graphically represents the comparison:
<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/Met_score_comparison.png" alt="img_1" class="inline"/>

#### Use the Jaccard index to quantify the modified cytosines

We can use Jaccard similarity coefficient to score the similarity between the regions that have been recognised by our tools as being modified. The Jaccard score measures the intersection over the union minus the intersection. So that, the closest the index is to one, the more the 2 datasets are similar.  
To do so, we can use bedtools. First we have to convert the frequency files to BED files (which is the format required from `bedtools`).

```{}
for i in ${OUT}/Met_frequency_*.tsv; do awk -F "\t" '{ if (NR > 1 ) print $1,$2,$2,$(NF-1) }' $i | sort -V | sed 's/ /\t/g' > ${i:0:-3}bedgraph ; done
```
Then we can calculate the Jaccard score for every pair of methods:
```{}
bedtools jaccard -a ${OUT}/Met_frequency_deepsignal.bedgraph -b ${OUT}/res/bisulfite.ENCFF835NTC.example.tsv
bedtools jaccard -a ${OUT}/Met_frequency_nanopolish.bedgraph -b ${OUT}/res/bisulfite.ENCFF835NTC.example.tsv
bedtools jaccard -a ${OUT}/Met_frequency_deepsignal.bedgraph -b ${OUT}/Met_frequency_nanopolish.bed
```
the output is:
```{}
intersection	union-intersection	jaccard 	n_intersections
```
the index we want to use to compare the regions is the third value (i.e. `jaccard`)

#### Compare the values assigned to the modified or not modified cytosines

We can also plot the similarity of the regions that have been recognised to be modified versus the one that have not been recognised. You can use the command:
```{}
${IN}/scripts/compare_met_regions.R -d ${OUT}/Met_frequency_deepsignal.tsv -n ${OUT}/Met_frequency_nanopolish.tsv -b ${IN}/res/bisulfite.ENCFF835NTC.example.tsv -o ${OUT}/methylation_region_comparison
```
Also this command print a table and a PDF file. The table report if a region has been recognised (TRUE) or not (FALSE) from a method/technique. It should look like: 
```{}
reference_regions	  deep_overlap	nano_overlap	bisul_overlap
chr20:5092445	      FALSE	        FALSE	        TRUE
chr20:5092454	      TRUE	        FALSE	        TRUE
chr20:5092455	      TRUE	        FALSE	        TRUE
chr20:5092482	      TRUE	        TRUE	        TRUE
chr20:5092483	      TRUE	        FALSE	        TRUE
chr20:5092485	      TRUE	        FALSE	        TRUE
chr20:5092486	      TRUE	        FALSE	        TRUE
```
and the PDF:
<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/met_region_comparison.png" alt="img_1" class="inline"/>

## Methylation - Visualisation {#Methylation-Visualization}

We could use IGV to visualise the bedgraph files we created before. Bedgraph files are essentially BED files with a 4th column that stores the value of that region. In the file we created the 4th column store the probabilty of the regions to be methylated. If we want to add the bisulphite as "reference methylation" we could convert it to bedgraph:
```{}
awk -F "\t" '{print $1,$2,$3,$11/100*$10}' ${IN}/res/bisulfite.ENCFF835NTC.example.tsv > ${OUT}/bisulfite.ENCFF835NTC.example.bedgraph
```
And then open all the files on IGV. It should looks something like:
<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/met_comparison_igv.png" alt="img_1" class="inline"/>
Do you think IGV tracks reflect the comparison plots we did before?


Alternatively, we could write a script to visualize the data. Scripts are a useful way to have consistency in your results and pictures. We can use an R script to plot the methylation peaks, called using the three different methods we saw during the course.

Open the script `${IN}/scripts/plot_methylation_gviz.R` and plot the methylation peaks using R. You should get something like:
<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/gviz_met.png" alt="img_1" class="inline"/>

## References
[Peng Ni, Neng Huang, Zhi Zhang, De-Peng Wang, Fan Liang, Yu Miao, Chuan-Le Xiao, Feng Luo, Jianxin Wang, DeepSignal: detecting DNA methylation state from Nanopore sequencing reads using deep-learning, Bioinformatics, Volume 35, Issue 22, 15 November 2019, Pages 4586–4595](https://doi.org/10.1093/bioinformatics/btz276)
[Simpson, Jared T., et al. “Detecting DNA cytosine methylation using nanopore sequencing.” nature methods 14.4 (2017): 407-410](https://doi.org/10.1038/nmeth.4184)
[The data for this practical have been taken from the Nanopolish GitHub page](https://github.com/jts/nanopolish)
