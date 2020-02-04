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
└── data
    └── 20180521_FAH88251
        └── reads
            ├── 0
            └── 1
```

This is what a raw directory structure will look like if you **dont** let the ONT software cal it for you!

Lets take a closer look at the reads/ directory. This is where the sequencer writes raw data.

```
ls ./data/20180521_FAH88251/reads/
```

You should see two directories: 0/ and 1/

ON machines write reads in batches, here batch size is 8000. This is to keep the number of files in each directory resonable for the computers filesystem

Take a look into a read directory!

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
h5ls data/20180521_FAH88251/reads/0 GXB01206_20180518_FAH88251_GA40000_mux_scan_MK_43023_read_103_ch_219_strand.fast5
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
/Raw/Reads/Read_103      Group
/Raw/Reads/Read_103/Signal Dataset {8483/Inf}
/UniqueGlobalKey         Group
/UniqueGlobalKey/channel_id Group
/UniqueGlobalKey/context_tags Group
/UniqueGlobalKey/tracking_id Group
```

We can also dump and view the entire contents of a FAST5 to text:

```
h5dump reads/0/GXB01206_20180518_FAH88225_GA50000_sequencing_run_CD3_92236_read_9998_ch_295_strand.fast5 | less
```

(**Hint 1**: use the mouse wheel to scroll up and down the file)

(**Hint 2**: press q to exit less, a text reading program)
