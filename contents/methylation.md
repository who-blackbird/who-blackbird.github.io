#  Oxford Nanopore - day 2 - Methylation

#### [Methylation - DeepSignal](#Methylation-DeepSignal)
#### [Methylation - Nanopolish](#Methylation-Nanopolish)
#### [Methylation - Comparison](#Methylation-Comparison) 
#### [Methylation - Visualisation](#Methylation-Visualization)

***

N.B.: Note these scripts use `course/` as working directory

## Methylation - DeepSignal {#Methylation-DeepSignal}

Deepsignal [(doi:10.1093/bioinformatics/btz276)](https://doi.org/10.1093/bioinformatics/btz276) is a Machine learning (ML) approach to identify changes in the resistance measured by the ONT that may be identified as a 5-MethylCytosin (5mC). As all ML approaches, DeepSignal requires a model. The model has been trained by the authors and you can download it with the software. DeepSignal takes ~ 24 hours to run on 4 GPU for a human genome 30x coverage.

Deepsignal is a quite straightforward tool. However, it needs:
+ Fast5 files to work (they contain the resistance information required to call 5mC).
+ if fast5 files have a Multi-fast5 format, they need to be converted to single fast5 (you can use the command [```multi_to_single_fast5```](https://github.com/nanoporetech/ont_fast5_api)).
+ It needs the ML model. One model has been trained by the Author. Alternatively, one can train a new model, especially if one has to investigate DNA/RNA modifications that are different from 5mC. To train a new model, one can use the command ```deepsignal train```. In this case, one needs a _training set_ and a _validation set_ (N.B. They need to be 2 different data sets).

#### Call 5mC modification

```{}
deepsignal call_mods --input_path methylation/res/fast5_files/ --model_path methylation/res/model.CpG.R9.4_1D.human_hx1.bn17.sn360/bn_17.sn_360.epoch_7.ckpt --result_file methylation/Met_deepsignal.tsv --reference_path methylation/res/reference.fasta --corrected_group methylation/RawGenomeCorrected_000 --nproc 10 --is_gpu no
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
python3 methylation/scripts/call_modification_frequency.py --input_path methylation/Met_deepsignal.tsv --result_file methylation/Met_frequency_deepsignal.tsv --prob_cf 0.6
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
nanopolish index -d methylation/res/fast5_files/ methylation/res/alignment_output.fastq
```
We need to align the reads to the reference genome. We can use mimimap2 for this and pipe the results to samtools to sort them. Ultimately, we need to index the bam file we created:
#### Align and sort fastq files
```{}
# Align the reads the the GRCh38 subset
minimap2 -a -x map-ont methylation/res/reference.fasta methylation/res/alignment_output.fastq | samtools sort -T tmp -o methylation/res/alignment_output_sorted.bam

# Index the sorted output
samtools index methylation/res/alignment_output_sorted.bam
```
Once done this preprocessing, we have all the required files to start the 5mC modification call. We use the command `call-methylation` to annotate the modification.
#### Call 5mC modification
```{}
nanopolish call-methylation -t 8 -r methylation/res/alignment_output.fastq -b methylation/res/alignment_output_sorted.bam -g methylation/res/reference.fasta -w "chr20:5,000,000-10,000,000" > methylation/Met_nanopolish.tsv
```
__NB.: This command will take about 10 minutes, so it is a good time to have a break and stretch a bit your.__

Tip:
`-t` specify the number of threads to use for this command. This value can be optimised on the number of cores in your CPU. To check how many threads you have, on linux, you can do:
```{}
lscpu | egrep 'Model name|Socket|Thread|NUMA|CPU\(s\)'
```
This command will see CPUs and threads which are the information you need, there is not an optimum number of threads, it really depends on the process you are running and it always need a bit of optimization.


#### Call 5mC modification frequency

```{}
methylation/scripts/nanopolish/scripts/calculate_methylation_frequency.py methylation/Met_nanopolish.tsv > methylation/Met_frequency_nanopolish.tsv
```

## Methylation - Comparison {#Methylation-Comparison}

#### Compare the cytosines that have been identified to be modified or not modified

```{}
methylation/scripts/compare_met_score.R -d methylation/Met_deepsignal.tsv -n methylation/Met_frequency_nanopolish.tsv -b methylation/res/bisulfite.ENCFF835NTC.example.tsv -o methylation/methylation_score_comparison
```

#### Use the Jaccard index to quantify the modified cytosines

We are using bedtools to calculate the Jaccard score.

Convert every frequency file to a BED file (required for bedtools)
```{}
for i in methylation/Met_frequency_*.tsv; do awk --field-separator="\t" '{ if (NR > 1 ) print $1,$2,$2,$(NF-1) }' $i | sort -V | sed 's/ /\t/g' > ${i:1:-3}bedgraph ; done
```

Calculate the Jaccard score
```{}
bedtools jaccard -a methylation/Met_frequency_deepsignal.bed -b methylation/res/bisulfite.ENCFF835NTC.example.tsv
bedtools jaccard -a methylation/Met_frequency_nanopolish.bed -b methylation/res/bisulfite.ENCFF835NTC.example.tsv
bedtools jaccard -a methylation/Met_frequency_deepsignal.bed -b methylation/Met_frequency_nanopolish.bed

```

#### Compare the values assigned to the modified or not modified cytosines

```{}
methylation/scripts/compare_met_regions.R -d methylation/Met_frequency_deepsignal.tsv -n methylation/Met_frequency_nanopolish.tsv -b methylation/res/bisulfite.ENCFF835NTC.example.tsv -o methylation/methylation_region_comparison
```

## Methylation - Visualisation {#Methylation-Visualization}

```{}

```
convert the `.Met_deepsignal.tsv` file to a BED file with this script and visualise it with the others tracks we have loaded earlier.
```{}

```

Alternatively, we could write a script to visualize the data. Scripts are a useful way to have consistency in your results and pictures. Here, we use an R script to plot the methylation peaks, called using the three different methods we saw during the course.

We load the packages we will need
```{}
library(Gviz)
library(tidyverse)
library(GenomicRanges)
```

We then define all the variables we are going to use:
```{}
# Define the files we need and their paths
deep <- "Desktop/NANOPORE_HPC/course/Met_frequency_deepsignal.bed"
nano <- "Desktop/NANOPORE_HPC/course/Met_frequency_nanopolish.bedgraph"
met <- "Desktop/NANOPORE_HPC/course/res/bisulfite.ENCFF835NTC.example.bed"
# Set the chromosomal coordinates to plot
chr <- "chr20"
start_reg <- 5000000
end_reg <- 6000000
genome_plot="hg38"
# Set the window number
window_plot = 250 # the final plot will display the regions divided in the number or windows we define here
```
Import the data:
```{} 
# Load the files as dataframe
deep_df <- read_delim(deep, delim = "\t", col_names = FALSE) %>%  
  filter(X2 < end_reg) # we filter out the regions that we don't want to plot (i.e. the ones starting  downstream our plot end) so that we work with a smaller dataset
nano_df <- read_delim(nano, delim = "\t", col_names = FALSE) %>%  
  filter(X2 < end_reg)
met_df <- read_delim(met, delim = "\t", col_names = FALSE) %>% 
  mutate(value = (X11 / 100 * X10) / max(X11 / 100 * X10) ) %>% # it calculates the methylation frequency and normalise it to 1
  filter(X2 < end_reg)
```
We want to see how the methylation correlates to the genes, possibly to their promoters. To download these info, it is convient to write down a few lines of code and get this info from ENSEMBL. If we are going to use this info quite often, it is convenient to save it and not download it every time, so that next time we can import it without downloading.
```{}
# Download from ENSMBL the gene symbols and coordinates
mart <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

attributes_to_extract <- c("chromosome_name", "start_position", "end_position", "external_gene_name")

values_to_filter <- paste(chr, format(start_reg, scientific = FALSE), format(end_reg, scientific = FALSE), sep = ":") %>% 
                    str_replace("chr","") # remove the 'chr' at the beginning. ENSEMBL doesn't use it 

gene_name <- biomaRt::getBM(attributes = attributes_to_extract, # Refer to the values you want to have in the final df
                                    filter="chromosomal_region", # Define what filter to apply. 
                                    values = values_to_filter, 
                                    mart = mart #Mart object to use
)  

gene_name$chromosome_name <-   str_replace(gene_name$chromosome_name, "^", "chr") #Put 'chr' on the chromosome

gene_name <- gene_name %>% 
  filter(!str_detect(external_gene_name, "LIN|\\.|-|SLC23A2")) # Remove long non coding RNA, MT genes and antisense RNA

# Transform the ENSEMBL data to a GenomicRange
gene_name_gr <- GenomicRanges::GRanges(seqnames = gene_name$chromosome_name, 
                                       ranges = IRanges(start = gene_name$start_position, 
                                                        end = gene_name$end_position) 
)
mcols(gene_name_gr)$gene_symbol <- gene_name$external_gene_name  
```
Ultimately, we can start to make the track we will use to plot.
```{}
# This line creates the chromosome(s) we want to plot, with giemsa banding and higligths the region of the cromosome we are plotting with a red box
ideoTrack <- IdeogramTrack(genome=genome_plot, chromosome=chr, start = start_reg, end = end_reg)
# This track gives the genomic coordinates we are plotting in Mb. a sort of ruler
gtrack <- GenomeAxisTrack()
# These 3 lines plot the methylation profile for deepsignal, nanopolish and bisulphite methods respectively.
plotTrack_deep <- DataTrack(range=deep_gr, genome=genome_plot, chromosome=chr, name="Methylation deepsignal", start = start_reg, end = end_reg, window = window_plot, type = c("a","hist"))
plotTrack_nano <- DataTrack(range=nano_gr, genome=genome_plot, chromosome=chr, name="Methylation nanopolish", start = start_reg, end = end_reg, window = window_plot, type = c("a","hist"))
plotTrack_met <- DataTrack(range=met_gr, genome=genome_plot, chromosome=chr, name="Methylation bisulphite", start = start_reg, end = end_reg, window = window_plot, type = c("a","hist"))

# This last line of script creates the graphic object with the plot we created so far.
#### plot the 
plotTracks(list(ideoTrack, gtrack, plotTrack_deep, plotTrack_nano, plotTrack_met, gene_track), chromosome=chr, start = start_reg, end = end_reg, background.title="white", col = "black", fontcolor = " black", cex.title = 1, col.title="black", col.border.title = "black", col.id="black", cex.axis=0.6, col.axis="black")

# to see all the graphic parameters go to https://rdrr.io/bioc/Gviz/man/settings.html
```
You should get something like:
<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/gviz_met.png" alt="img_1" class="inline"/>