# ABySS-Connector: Connecting Paired Sequences with a Bloom Filter de Bruijn Graph
#### _Note: Submission Requirements for HiTSeq 2014_

Full paper submissions to HiTSeq 2014 have the same requirements as "Original Paper" submissions to "Bioinformatics", as described at http://www.oxfordjournals.org/our_journals/bioinformatics/for_authors/general.html. Briefly, the requirements are:

* Submission deadline: April 6, 2014
* Preferred submission formats: Word 2003 for document, .TIF or .EPS for images
* Up to 7 pages in length (~ 5000 words)
* The research must be demonstrated with real biological data, not just simulated data
* Required sections:
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


