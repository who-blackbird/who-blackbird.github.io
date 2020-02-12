# Quality control

**The run completed successfully - so now what?**

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/qc.starting_point.png" alt="img_1" class="inline"/>

### In this section we will cover:

1. Signal to sequence workflow
2. Run directory structure
3. FAST5 files
4. Basecalling
5. Basecalling results

#### You will learn to:

7. Take a run summary file and extract QC data
8. Interpret this data

## 1. Singal to sequence workflow

The below image represents the process of translating raw electrical signal data from an ON sequencer to DNA sequence.

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/qc.basecalling.png" alt="img_2" class="inline"/>

We are going to learn how to do this!

But... before we get too technical, lets take a look at how Nanopore's machines write data!

## 2. Run directory structure

All nanopore data is written to a specific run directory. Open up a terminal and change to the example directory below…

```
cd /home/participant/Course_Materials/data/quality_control
```

... and check it out using tree!

```
tree -d .
```

_(**Hint**: "**.**" is a shortcut for the current directory you are in)_  
_(**Hint 2**: "**-d**" means the tree command lists directories only)_  
_(**Hint 3**: at home mac users might need to install the tree command!)_

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

Lets take a closer look at the **reads/** directory. This is where the sequencer writes raw data:

```
ls 20180830_PAD01151/data/reads/
```

You should see two directories: **0/** and **1/**.

ON machines write reads in batches, here batch size is **100** files. Normally, batch size would be around **5000** to **8000**. This is to keep the number of files in each directory resonable for the computers filesystem.

Take a look into read directory **0/**:

```
ls ./data/20180521_FAH88251/reads/0/
```

You should see a list of files ending with **.fast5**.

## 3. FAST5 files

Nanopore writes read data to a file format they call **FAST5**

1 read = 1 **.fast5** file

**FAST5** files are infact **HDF5** files. These are compressed binary files which store data in a structured way, allowing fast random access to subsets of the data.

This is where electronic signal data from the sequencer is storred.

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/qc.squiggle_data.png" alt="img_2" class="inline"/>

Lets look at the structure of a **FAST5** file using the **h5ls** command:

```
h5ls 20180830_PAD01151/data/reads/0/PCT0045_20180830_0004A30B002402F4_2_E3_H3_sequencing_run_Cam_6_90217_read_100_ch_1023_strand.fast5
```

This shows the top level data keys:

```
PreviousReadInfo         Group
Raw                      Group
UniqueGlobalKey          Group
```

We can view the subkeys by recursivley (**-r**) listing the file:

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

We can also dump and view the entire contents of a **.fast5** to text:

```
h5dump 20180830_PAD01151/data/reads/0/PCT0045_20180830_0004A30B002402F4_2_E3_H3_sequencing_run_Cam_6_90217_read_100_ch_1023_strand.fast5 | less
```

(**Hint 1**: use the up and down arrows to scroll through the file)  
(**Hint 2**: press **q** to exit less, a text reading program)

## 4. Basecalling

This is the process of translating raw electrical signal data from an ON sequencer to DNA sequence. Basecalling is a critical step in the analysis workflow as poor basecalling results in poor sequence data.

Many basecallers exist - but for now we will be using **Guppy v3.4.5** developed by Oxford Nanopore.

For basecalling it is important to know which **Flow Cell** and **Library Prep Kit** was used.

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

**Interesting fact**: You have just started up a machine learning algorithm. Guppy, alongside almost all current nanopore basecallers have a neural network at their core.

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

## 5. Basecalling results

Lets take a look inside the new **basecalled/** directory we just created:

```
ls -l 20180830_PAD01151/data/basecalled

```

```
total 4024
-rw-r----- 1 nickgleadall staff 1797686 4 Feb 15:52 fastq_runid_5f42e62fe05de2c3f4e892def6936adeb53cf710_0_0.fastq
-rw-r----- 1 nickgleadall staff 121579 4 Feb 15:52 guppy_basecaller_log-2020-02-04_15-47-12.log
-rw-r----- 1 nickgleadall staff 57128 4 Feb 15:52 sequencing_summary.txt
-rw-r----- 1 nickgleadall staff 65669 4 Feb 15:52 sequencing_telemetry.js
drwxr-x--- 4 nickgleadall staff 128 4 Feb 15:47 workspace
```

Our freshly called "sequence" data now sits inside **fastq_runid_5f42e62fe05de2c3f4e892def6936adeb53cf710_0_0.fastq**.

Lets take a quick look at it:

```
less 20180830_PAD01151/data/basecalled/fastq_runid_5f42e62fe05de2c3f4e892def6936adeb53cf710_0_0.fastq
```

This is a standard format fastq file which can now be analysed downstream. We will cover fastq's in greater detail later!  
(_**Hint**: press q to exit less, that text reading program from before!_)

Another important file is **sequencing_summary.txt**. This contains **A LOT** of data about the sequencing run. Lets take a look.

```
less 20180830_PAD01151/data/basecalled/sequencing_summary.txt
```

## Quick breather!!!

What **should** I know at this point?:

- We are now at the point where we understand how raw squiggle data is converted into to A's T's C's and G's.

- We can take some .fast5 files, run a basecaller, and find the results.

- We have also learned that there are many different tools for basecalling the data. However, for general use Guppy from Oxford Nanopore is fine.

- It is also recommended that you let the sequencing software handle this automatically for you.

What **dont** we know:

- How can we access the data contained in **sequencing_summary.txt**.

- How do we interpret this data?

Now, back to work!

## 6. Take a run summary file and extract QC data

Many tools for run QC have been developed by the Nanopore sequening community. We will now learn how to use one of them!

Before we begin, lets change directory.

```
cd ~/Desktop/quality_control/QC_script
```

We can then open up a file browser to make things easy.

```
open .
```

This should open up a directory, and we can see that there are multiple files in this "package".

The first file we are interested in is the **config.yaml** file.
_**Hint**: dont worry .yaml is just a fancy text file!_

This file contains the configuration for the QC plotting tool.

On opening the file you should see the following configuration lines at the top:

```
inputFile: "./RawData/lambda_sequencing_summary.txt.bz2"
barcodeFile: "./RawData/lambda_barcoding_summary.txt.bz2"
basecaller: "Guppy 2.1.3"
flowcellId: "EXAMPLE_FLOWCELL"
tutorialText: FALSE
```

Notice that the input file points to a version of the **sequencing_summary.txt**

We are going to leave this alone for now!

The second file we are interested in is the **Nanopore_SumStatQC.Rmd** file. Click it and Rstudio should open up.

Press the "Knitr" button at the top of the page and wait for magic!

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/qc.rhelp.png" alt="img_2" class="inline"/>

_**Disclaimer**: If magic does not happen, R should have printed a .html file - open it by pasting the below in the terminal!_

```
open Nanopore_SumStatQC.html
```

## 7. Interpret this data

We can read through the report to get an understanding of QC metrics!

## 8. Bonus step - a real run!

Lambda controls are so boring. Lets take a look at the QC report of a PromethION run from a patient with a rare bleeding disorder!!!

Tree the following directory:

```
tree ../example_runs/2.promethion_run
```

_Hint: "**..** is a shortcut for "one directory up"_

You should see:

```
../example_runs/2.promethion_run
└── VWD1108
    └── 20191202_1317_1-E11-H11_PAE13924_73938d1f
        ├── PCT0099_20191202_131743_PAE13924_promethion_sequencing_run_VWD1108_duty_time.csv
        ├── PCT0099_20191202_131743_PAE13924_promethion_sequencing_run_VWD1108_sequencing_summary.txt
        ├── PCT0099_20191202_131743_PAE13924_promethion_sequencing_run_VWD1108_throughput.csv
        ├── fast5_fail
        │   ├── PAE13924_3a68f5baa685d37c37fb291fde3b6f0c6b02573c_29.fast5
        │   └── PAE13924_3a68f5baa685d37c37fb291fde3b6f0c6b02573c_86.fast5
        ├── fast5_pass
        │   ├── PAE13924_3a68f5baa685d37c37fb291fde3b6f0c6b02573c_787.fast5
        │   └── PAE13924_3a68f5baa685d37c37fb291fde3b6f0c6b02573c_991.fast5
        ├── fast5_skip
        │   └── PAE13924_3a68f5baa685d37c37fb291fde3b6f0c6b02573c_0.fast5
        ├── fastq_fail
        │   ├── PAE13924_3a68f5baa685d37c37fb291fde3b6f0c6b02573c_433.fastq
        │   └── PAE13924_3a68f5baa685d37c37fb291fde3b6f0c6b02573c_490.fastq
        ├── fastq_pass
        │   ├── PAE13924_3a68f5baa685d37c37fb291fde3b6f0c6b02573c_816.fastq
        │   └── PAE13924_3a68f5baa685d37c37fb291fde3b6f0c6b02573c_991.fastq
        ├── final_summary.txt
        ├── report.md
        └── report.pdf
```

This is what a "run" directory looks like when you let the sequencer take care of data processing.

You will notice that **.fast5** and **.fastq** files have beeen automatically separated into **pass** and **fail** directories.

There is also a "report.pdf" which the sequencer prints off, however the R script we used before gives us much more information (feel free to take a look though!).

Lets produce the report for this run now. Open up the **config.yaml** file and edit the lines:

```
inputFile: "./RawData/lambda_sequencing_summary.txt.bz2"
barcodeFile: "./RawData/lambda_barcoding_summary.txt.bz2"
basecaller: "Guppy 2.1.3"
flowcellId: "Lambda test"
tutorialText: FALSE
```

So that they read:

```
inputFile: "../example_runs/2.promethion_run/VWD1108/20191202_1317_1-E11-H11_PAE13924_73938d1f/PCT0099_20191202_131743_PAE13924_promethion_sequencing_run_VWD1108_sequencing_summary.txt"
barcodeFile: ""
basecaller: "Guppy 2.1.3"
flowcellId: "PAE24233"
tutorialText: FALSE
```

_**Hint**: copy and paste._

Save the file.

Now press the "**Knitr**" button again the top of the page and wait for some _...slightly slow..._ magic!

You can print these reports to **.pdf** used the standard browser method if required.

If you want to make youre own QC pipeline or just see how this one works, read through the **Nanopore_SumStatQC.Rmd** script.

## Final Comment

This was a VERY basic overview of nanopore data analysis. Below is a diagram showing the parts of an overall workflow we have covered so far.

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/qc.so_far.png" alt="img_2" class="inline"/>

## References & Reading

Excellent comparison of basecallers: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1727-y

Oxford Nanopore downloads: https://community.nanoporetech.com/downloads

Oxford Nanopore tutorials: https://community.nanoporetech.com/knowledge/bioinformatics/

Specific QC tutorial: https://github.com/nanoporetech/ont_tutorial_basicqc

How are Phred scores calculated: https://en.wikipedia.org/wiki/Phred_quality_score
