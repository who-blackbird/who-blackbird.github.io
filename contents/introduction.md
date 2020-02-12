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
  <img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/intro.genome_sizes.png" alt="img_2" class="inline"/>
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
  <img src="//raw.githubusercontent.com/who-blackbird/who-blackbird.github.io/master/images/intro.lrs.png" alt="img_2" class="inline"/>
</p>

The two best known producers of ‘true’ long-read sequencing technologies are Pacific Biosciences (**PacBio**) and Oxford Nanopore Technologies (**Nanopore**). Bother offer platforms that directly sequence single molecules of **DNA or RNA** in '**real-time**' and which can produce significantly longer reads than current short-read platofrms.

## 2. SMRT sequencing

## References
