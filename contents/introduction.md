# Introduction

### In this section we will cover:

1. What is long range sequencing?
2. SMRT sequencing
3. Nanopore sequencing
4. Different ONT sequencers
5. Advantages of this method
6. Disadvantages of this method
7. Why have we been using nanopore?
8. A standard workflow

## 1. What is long read sequencing?

The genomes of most organisms are too big to be sequenced in one long continuous string.

<p align="center">
  <img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/intro.genome_sizes.png" alt="img_1" class="inline"/>
</p>

In order to sequence these massive genomes, **short-read sequencing** was developed. Here DNA is broken into millions of tiny fragments that are amplified (copied) and then sequenced.

This process produces millions of short sequences or "**reads**" which are typically **75-300bp** in length. Reads are then pieced together computationally.

<ins>**Short-read workflow**</ins>

<p align="center">
  <img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/intro.srs.png" alt="img_2" class="inline"/>
</p>

In contrast, **long-read sequencing** allows us to retrive much longer sequence reads, typically **>10,000bp** in length. In fact, some long read sequenceing systems have been producing reads over **2,000,000bp** (2 megabases) in length.

Additionally long-read methodologies typically sequence single molecules of DNA directly and in real-time, foregoing the need for amplification.

<ins>**Long-read workflow**</ins>

<p align="center">
  <img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/intro.lrs.png" alt="img_3" class="inline"/>
</p>

The two best known producers of ‘true’ long-read sequencing technologies are Pacific Biosciences (**PacBio**) and Oxford Nanopore Technologies (**Nanopore**). Bother offer platforms that directly sequence single molecules of **DNA or RNA** in '**real-time**' and which can produce significantly longer reads than current short-read platofrms.

## 2. SMRT sequencing

Single Molecule, Real Time (**SMRT**) sequencing the name of **PacBio**'s methodology.

SMRT sequencing workflow:

1. DNA is extracted from the sample using standard methods
2. DNA is then fragmented
3. Sequencing adapters are ligated to the fragments creating circular templates.

<p align="center">
  <img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/intro.circular.png" alt="img_4" class="inline"/>
</p>

4.  Primer and polymerase are then added to the sequencing library.
5.  The library is then loaded into a "**SMRT cell**" containing millions of 'wells' called Zero Mode Wave-guides (**ZMW**'s).
6.  A single molecule of DNA is immobilised inside a ZMW.
7.  Sequencing reaction now takes place. As the DNA polymerase incorporates nuclotides light is emitted and measured.

<p align="center">
  <img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/intro.smrt.jpg" alt="img_5" class="inline"/>
</p>

**PacBio** **SMRT** sequencing has two modes:

- **Circular Consensus Sequencing (CCS)**: Sequencies a smaller circular template over and over to produce a high quality read.
- **Continuous Long Read Sequencing (CLR)**: Produces the longest read by sequencing long templates.

Original Publication: Eid, J., et al. (2009). [Real-time DNA sequencing from single polymerase molecules](http://dx.doi.org/10.1126/science.1162986). Science, 323(5910), 133–138.

------ **Disclaimer** ------:

We have never actually used or worked with **PacBio** sequencers/data. From now on this course will focus on **Nanopore** sequencing!

If you would like more information about PacBio, I suggest you contact Professor Steven Marsh or Dr. James Robinson who are based at the [Anthony Nolan Research Institute in London](https://www.anthonynolan.org/clinicians-and-researchers/anthony-nolan-research-institute/hla-informatics-group).

They have been using **PacBio** sequencing for HLA typing over the past few years and have expert knowledge on this system.

## References
