# 3. Quality control

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/qc.starting_point.png" alt="img_1" class="inline"/>

**The run completed successfully - so now what?**

### In this section we will cover:

1. Signal to sequence workflow
2. Run directory structure
3. FAST5 files
4. Basecalling
5. Fastq

#### You will learn to:

7. Take a run summary file and extract QC data
8. Interpret this data

## 1. Singal to sequence workflow

The below image represents the process of translating raw electrical signal data from an ON sequencer to DNA sequence.

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/qc.basecalling.png" alt="img_2" class="inline"/>

We are going to learn how to do this!!

But... before we get too technical, lets take a look at how Nanopore's machines write data!

## 2. Run directory structure

All nanopore data is written to a specific run directory. Open up a terminal and change to the example directory below…

```
cd /Users/nickgleadall/Desktop/quality_control/example_runs/1.raw_data
```

... and check it out using tree!

```
tree -d .
```

_note: "**.**" is a shortcut for the current directory you are in_

_note2: "**-d**" means the tree command lists directories only_

_note3: at home mac users might need to install the tree command!_

You should see something like this:

```
.
└── 20180830_PAD01151
    └── data
        └── reads
            ├── 0
            └── 1

5 directories
```

This is what a raw directory structure will look like if you **dont** let the ONT software cal it for you!

Lets take a closer look at the reads/ directory. This is where the sequencer writes raw data.

```
ls 20180830_PAD01151/data/reads/
```

You should see two directories: 0/ and 1/

ON machines write reads in batches, here batch size is 100 files. Normally, batch size would be around 5000 to 8000. This is to keep the number of files in each directory resonable for the computers filesystem.

Take a look into read directory 0/!

```
ls ./data/20180521_FAH88251/reads/0/
```

You should see a **HUGE** list of files ending with .fast5

## 3. FAST5 files

Nanopore writes read data to a file format they call FAST5

1 read = 1 .fast5 file

FAST5 files are infact HDF5 files. These are compressed binary files which store data in a structured way, allowing fast random access to subsets of the data.

This is where electronic signal data from the sequencer is storred.

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/qc.basecalling.png" alt="img_2" class="inline"/>

Lets look at the structure of a FAST5 file using h5ls:

```
h5ls 20180830_PAD01151/data/reads/0/PCT0045_20180830_0004A30B002402F4_2_E3_H3_sequencing_run_Cam_6_90217_read_100_ch_1023_strand.fast5
```

This shows the top level data keys:

```
PreviousReadInfo         Group
Raw                      Group
UniqueGlobalKey          Group
```

We can view the subkeys by recursivley (-r) listing the file:

```
h5ls -r data/20180521_FAH88251/reads/0 GXB01206_20180518_FAH88251_GA40000_mux_scan_MK_43023_read_103_ch_219_strand.fast5
```

Outputs:

```
/                        Group
/PreviousReadInfo        Group
/Raw                     Group
/Raw/Reads               Group
/Raw/Reads/Read_100      Group
/Raw/Reads/Read_100/Signal Dataset {25746/Inf}
/UniqueGlobalKey         Group
/UniqueGlobalKey/channel_id Group
/UniqueGlobalKey/context_tags Group
/UniqueGlobalKey/tracking_id Group
```

We can also dump and view the entire contents of a FAST5 to text:

```
h5dump 20180830_PAD01151/data/reads/0/PCT0045_20180830_0004A30B002402F4_2_E3_H3_sequencing_run_Cam_6_90217_read_100_ch_1023_strand.fast5 | less
```

(**Hint 1**: use the mouse wheel to scroll up and down the file)

(**Hint 2**: press q to exit less, a text reading program)

## 4. Basecalling

This is the process of translating raw electrical signal data from an ON sequencer to DNA sequence. Basecalling is a critical step in the analysis workflow as poor basecalling results in poor sequence data.

Many basecallers exist - but for now we will be using Guppy v3.4.5developed by ON.

For basecalling it is important to know which Flow Cell and Library Prep Kit was used.

To see all the combinations which Guppy can handle use:

```
guppy_basecaller --print_workflows | head -n 10
```

This will list all the supported combinations of **flowcell** and **library prep kit** that Guppy can basecall.

```
Available flowcell + kit combinations are:
flowcell   kit        barcoding config_name
FLO-FLG001 SQK-RNA003           rna_r9.4.1_120bps_hac
FLO-MIN106 SQK-RNA003           rna_r9.4.1_120bps_hac
FLO-MIN107 SQK-RNA003           rna_r9.4.1_120bps_hac
FLO-MIN107 SQK-DCS108           dna_r9.5_450bps
FLO-MIN107 SQK-DCS109           dna_r9.5_450bps
FLO-MIN107 SQK-LRK001           dna_r9.5_450bps
FLO-MIN107 SQK-LSK108           dna_r9.5_450bps
FLO-MIN107 SQK-LSK109           dna_r9.5_450bps
```

So lets basecall!

```
guppy_basecaller --input_path 20180830_PAD01151/data/reads/ --recursive -s 20180830_PAD01151/data/basecalled --fast5_out --flowcell FLO-PRO002 --kit SQK-LSK109
```

(**Interesting fact**: You have just started up a machine learning algorithm. Guppy, alongside almost all current nanopore basecallers have a neural network at their core)

## 4.i. A comment about basecallers

As previously mentioned many basecallers are available.

The main performance marker of a basecaller that we care about is the overall **Assembly Identity** or how much a final alignment matches the reference.

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/qc.basecaller_comparison.png" alt="img_2" class="inline"/>

We can also take a look at the **assembly length bias** which tells us if a given basecaller is prone to reference insertions or deletions.

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/qc.assembly_length_bias.png" alt="img_2" class="inline"/>

_Further reading_: A great comparison of basecallers exists [here](https://github.com/rrwick/Basecalling-comparison/)

## 4.ii. So what do we recommend?

**In our group we have never used anything other than the most up to date ON basecalling algorithm (i.e. Guppy).**

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/qc.workflow.png" alt="img_2" class="inline"/>
