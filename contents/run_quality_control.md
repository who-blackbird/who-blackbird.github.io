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
cd ~/Desktop/quality_control/example_runs/example_1_20191202_PAE13924
```

... and check it out using tree!

```
tree .
```

_note: **.** is a shortcut for the current directory you are in_
_note2: at home mac users might need to install the tree command!_

#########################################

This process is called **basecalling**.

Basecalling is a critical step in the analysis workflow as poor basecalling makes poor sequence data.

#### Before we proceed...

As this is a 'introduction' course, we recommened that you enable 'on board' or 'real time' basecalling. In laymans terms, this means letting nanopores software automatically base call for you.

Infact, this is what we use as a group.
