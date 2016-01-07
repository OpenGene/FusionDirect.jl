# FusionDirect

[![Build Status](https://travis-ci.org/OpenGene/FusionDirect.jl.svg?branch=master)](https://travis-ci.org/OpenGene/FusionDirect.jl)

# FusionDirect
detect gene fusion directly from fastq files

## Features
* no alignment needed, it just reads fastq files of pair sequencing
* output fusion pattern (gene and position), along with the reads support this fusion
* ultra sensitive, comparing to delly, factera or other tools
* output file is a standard fasta file, which can be used to verify fusions using blast or other tools
* very suitable for detecting fusions from cancer target sequencing data (exom seq or panel seq)

## Install
```julia
# from Julia REPL
Pkg.clone("https://github.com/OpenGene/FusionDirect.jl.git")
```

## Use as a package
```julia
using FusionDirect

# the reference folder, which contains chr1.fa, chr2fa...
# can be downloaded from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
ref = "/opt/ref/hg19chr"
# a gene list with their coordination intervals, see the example bed files in data folder
bed = Pkg.dir("FusionDirect") * "/data/test_panel.bed"
read1 = "R1.fq.gz"
read2 = "R2.fq.gz"
detect(ref, bed, read1, read2)
```

## Use FusionDirect as a standalone script from commandline
copy src/fusion.jl to anywhere you want, run
```shell
julia fusion.jl -f <REF> -b <BED> -l <READ1> -r <READ2> > output.fa
```

## Sample output
```
#Fusion:ALK-EML4
>ALK|-|chr2:30144237_EML4|+|chr2:42396619/1
CCACCCAAAGCCGCGGGCGCTGATGATGGGTGAGGAGGGGGCGGCAAGATTTCGGGCGCCCCTGCCCTGAACGCCCTCAGCTGCTGCCGCCGGGGCCGCTCCGGAGGCGGGAGCCGGTAGCCGAGCCGGGCGACCTAGAGAACGAGCG
>ALK|-|chr2:30144237_EML4|+|chr2:42396619/2
CACTTCATTCAGCGGACCGACAGAGTGGCCGACGCTGAGCCTGACCCGCTCGTTCTCTAGGTCGCCCGGCTCGGCTACCGGCTCCCGCCTCCGGAGCGGCCCCGGCGGCAGCAGCTGAGGGCGTTCAGGGCAGGGGCGCCCGAAATC
>EML4|-|chr2:42396619_ALK|+|chr2:30144237/1
CACTTCATTCAGCGGACCGACAGAGTGGCCGACGCTGAGCCTGACCCGCTCGTTCTCTAGGTCGCCCGGCTCGGCTACCGGCTCCCGCCTCCGGAGCGGCCCCGGCGGCAGCAGCTGAGGGCGTTCAGGGCAGGGGCGCCCGAAATC
>EML4|-|chr2:42396619_ALK|+|chr2:30144237/2
CCACCCAAAGCCGCGGGCGCTGATGATGGGTGAGGAGGGGGCGGCAAGATTTCGGGCGCCCCTGCCCTGAACGCCCTCAGCTGCTGCCGCCGGGGCCGCTCCGGAGGCGGGAGCCGGTAGCCGAGCCGGGCGACCTAGAGAACGAGCG
```
