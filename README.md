# ABySS-Connector: Connecting Paired Sequences with a Bloom Filter de Bruijn Graph
#### _Note: Submission Requirements for HiTSeq 2014_

Full paper submissions to HiTSeq 2014 have the same requirements as "Original Paper" submissions to "Bioinformatics", as described at http://www.oxfordjournals.org/our_journals/bioinformatics/for_authors/general.html. Briefly, the requirements are:

* Submission deadline: April 6, 2014
* Preferred submission formats: Word 2003 for document, .TIF or .EPS for images
* Up to 7 pages in length (~ 5000 words)
* The research must be demonstrated with real biological data, not just simulated data
* Recommended Structure:
  * Title page
  * Structured Abstract
    * Motivation
    * Results
    * Availability and Implementation
    * Contact
    * Supplementary Information
  * Introduction
  * System and methods
  * Algorithm
  * Implementation
  * Discussion
  * Acknowledgements
  * Funding
  * References

## Abstract

Paired-end sequencing yields a read from each end of a DNA molecule, typically leaving a gap of unsequenced nucleotides in the middle of the fragment. We have developed ABySS-Connector, a software tool to fill in the nucleotides of the unsequenced gap by navigating a de Bruijn Graph to find a path between the two reads and connect the pair. ABySS-Connector represents the de Bruijn graph using a Bloom filter, a probabilistic and memory-efficient data structure that represents a set. Our implementation is able to store the de Bruijn graph using a mean 1.5 bytes of memory per k-mer, a marked improvement over the typical hash table data structure. The memory usage per k-mer is independent of k, enabling its application to larger genomes. The use of a Bloom filter to represent a de Bruijn graph has previously been described for genome sequence assembly, a task which benefits from a second non-probabilistic data structure to enumerate the critical false positives. We observe that this additional data structure is unnecessary for connecting reads, reducing memory requirements. The de Bruijn graph of the 20-gigabase white spruce genome sequencing data, for example, can be represented in 40 gigabytes. Those k-mers observed only once are usually erroneous, and so discarded by using a counting Bloom filter. Constructing the Bloom filter is parallelized and distributed over multiple machines, and
connecting the reads is likewise parallelized and distributed. ABySS-Connector is expected to have broad applications in genomic analysis, including read alignment, sequence assembly and haplotype variant calling.

## Introduction
## Method

ABySS-Connector determines the unknown nucleotide sequence in the gap between two paired-end reads by finding a connecting path through a de Bruijn graph. If a suitable path or set of paths can be found, the sequence corresponding to the connecting path(s) is inserted between the two reads and the full sequence for the DNA fragment is output as a single "pseudoread".

ABySS-Connector operates in three steps.  In the first step, the kmers from all of the paired end reads are loaded into a bloom filter representing the de Bruijn graph.  In the second step, a bidirectional graph search is carried out between each pair of reads to find connecting paths.  In the third and final step, a consensus sequence is constructed for each read pair which has multiple connecting paths, or the connecting sequence is used unmodified in the case of a unique path.

### Loading the Bloom Filter de Bruijn Graph
### Finding Connecting Paths within the de Bruijn Graph
### Reconciling Alternate Paths between Read Pairs

## Results

### Synthetic Data: Human Chromosome 21
#### Coverage of Reference by Connected Read Pairs
#### Correctness of Generated Sequences

### Real Data: Illumina Sequencing for Human Individual NA19238

Data: https://www.ebi.ac.uk/ena/data/view/PRJEB4252

#### Coverage of Reference by Connected Read Pairs
#### Correctness of Generated Sequences


