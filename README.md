# pairsnp-c++

[![pairsnp-CI](https://github.com/gtonkinhill/pairsnp-cpp/workflows/pairsnp-cpp-CI/badge.svg)](https://github.com/gtonkinhill/pairsnp-cpp/actions)

## Installation

The c++ version can be installed manually or with conda as

```
conda install -c bioconda pairsnp
```


The c++ code can be built by running

```
cd ./pairsnp-cpp/
make
make install
```

If your compiler does not have support for open-mp you can run

```
make ompoff=1
```

## Quick Start

The c++ version can be run from the command line as

```
pairsnp -c msa.fasta > output.csv
```

additional options include

```
SYNOPSIS
  Fast pairwise SNP distance matrices
USAGE
  pairsnp [options] alignment.fasta[.gz] > matrix.csv
OPTIONS
  -h	Show this help
  -v	Print version and exit
  -s	Output in sparse matrix form (i,j,distance).
  -d	Distance threshold for sparse output. Only distances <= d will be returned.
  -k	Will on return the k nearest neighbours for each sample in sparse output.
  -c	Output CSV instead of TSV
  -i	Output sequence index inplace of header (useful for downstream processing)
  -t	Number of threads to use (default=1)
```
