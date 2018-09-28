# pairsnp-cpp


## Installation

The c++ version can be installed manually or with conda as

```
conda install -c gtonkinhill pairsnp
```


The c++ code relies on a recent version of Armadillo (currently tested on v8.6) and can be built by running

```
cd ./pairsnp/cpp/
./configure
make
make install
```

The majority of time is spend doing sparse matrix multiplications so linking to a parallelised library for this is likely to improve performance further.

At the moment there you may need to run `touch ./cpp/*` before compiling to avoid some issues with time stamps.


## Quick Start

The c++ version can be run from the command line as

```
pairsnp -c msa.fasta > output.csv
```

additional options include

```
SYNOPSIS
  Pairwise SNP similarity and distance matrices using fast matrix algerbra libraries
USAGE
  pairsnp [options] alignment.fasta[.gz] > matrix.csv
OPTIONS
  -h	Show this help
  -v	Print version and exit
  -s	Find the similarity matrix
  -c	Output CSV instead of TSV
  -n	Count comparisons with Ns (off by default)
  -b	Blank top left corner cell instead of 'pairsnp 0.0.1'
```