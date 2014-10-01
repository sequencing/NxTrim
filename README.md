nxtrim: Software to remove Nextera Mate Pair adapters and categorise reads according to the orientation implied by the adapter location.  This software is not commercially supported.

Copyright (c) 2014, Illumina, Inc. All rights reserved.

This software is provided under the terms and conditions of the BSD 2-Clause License

You should have received a copy of the BSD 2-Clause License along with this program. If not, see https://github.com/sequencing/licenses/.

####Dependencies

BOOST - we use Boost 1.55.0 but most recent versions are probably fine

####Installation
```
git clone git@github.com:sequencing/NxTrim.git
cd NxTrim
make
./nxtrim
```

You will also need to point the BOOST_ROOT environment variable at your boost installation eg.
``
export BOOST_ROOT=/your/boost/installation
``

####Usage
Trim the data:
```
nxtrim -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -O sample 
```

Assemble with [Velvet](https://www.ebi.ac.uk/~zerbino/velvet/):
```
velveth output_dir 55 -short -fastq.gz sample.se.fastq.gz -shortPaired2 -fastq.gz sample.pe.fastq.gz -shortPaired3 -fastq.gz sample.mp.fastq.gz -shortPaired4 -fastq.gz sample.unknown.fastq.gz
velvetg output_dir -exp_cov auto -cov_cutoff auto -shortMatePaired4 yes
```

Assemble with [SPAdes](http://bioinf.spbau.ru/spades):
```
cat sample.mp.fastq.gz sample.unknown.fastq.gz > sample.allmp.fastq.gz
spades.py -k 21,33,55,77 -t 4  --pe1-s sample.se.fastq.gz --pe2-12 sample.pe.fastq.gz --hqmp3-12 sample.allmp.fastq.gz --hqmp3-fr  -o output_dir
```
Note we concatenate the unknown/mp libraries for SPAdes.  This command is suitable for 2x151bp data, if you have 2x251bp then use `-k 21,33,55,77,127`.  

####Details:

The default behaviour expects raw fastq files from a Nextera Mate-Pair library kit in Reverse-Forward orientation.  Based on the location of the Nextera adapter sequence (if detected), nxtrim produces four different "virtual libraries":

* mp: read pairs that are large insert-size mate-pairs
* pe: read pairs that are short insert-sze paired-end reads
* se: single reads 
* unknown: a library of read-pairs that are mostly large-insert mate-pair, but possibly contain a small proportion of paired end contaminants

The trimmer will reverse-complement the reads such that the resulting libraries will be in Forward-Reverse orientation, this reverse-complementing can be disabled via the --norc flag.

####Example data:

https://basespace.illumina.com/s/TXv32Ve6wTl9

Free registration required.

####References:

http://res.illumina.com/documents/products/technotes/technote_nextera_matepair_data_processing.pdf

http://res.illumina.com/documents/products/appnotes/appnote-nextera-mate-pair-bacteria.pdf
