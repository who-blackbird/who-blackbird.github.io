# 3. Quality control

<img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/qc.starting_point.png" alt="img_1" class="inline"/>

## The run completed successfully - so now what?

### In this section we will cover:

1. Signal to sequence workflow
2. Run directory structure
3. Raw "fast5" data
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
cd /Users/nickgleadall/Desktop/quality_control/example_runs/example_run_1
```

... and check it out using tree!

```
tree .
```

_note: **.** is a shortcut for the current directory you are in_
_note2: at home mac users might need to install the tree command!_

You should see something like this:

```
.
└── example_sample
    └── 20191202_1317_1-E11-H11_PAE13924_73938d1f
        └── fast5
            ├── PAE13924_3a68f5baa685d37c37fb291fde3b6f0c6b02573c_0.fast5
            ├── PAE13924_3a68f5baa685d37c37fb291fde3b6f0c6b02573c_29.fast5
            ├── PAE13924_3a68f5baa685d37c37fb291fde3b6f0c6b02573c_787.fast5
            ├── PAE13924_3a68f5baa685d37c37fb291fde3b6f0c6b02573c_86.fast5
            └── PAE13924_3a68f5baa685d37c37fb291fde3b6f0c6b02573c_991.fast5

3 directories, 5 files
```

#########################################

This process is called **basecalling**.

Basecalling is a critical step in the analysis workflow as poor basecalling makes poor sequence data.

#### Before we proceed...

As this is a 'introduction' course, we recommened that you enable 'on board' or 'real time' basecalling. In laymans terms, this means letting nanopores software automatically base call for you.

Infact, this is what we use as a group.

## References & Links

Nanopore software: https://community.nanoporetech.com/downloads
