nxtrim: Software to remove Nextera Mate Pair adapters and categorise reads according to the orientation implied by the adapter location.  This software is not commercially supported.

Copyright (c) 2014, Illumina, Inc. All rights reserved.

This software is provided under the terms and conditions of the BSD 2-Clause License

You should have received a copy of the BSD 2-Clause License along with this program. If not, see https://github.com/sequencing/licenses/.

####Dependencies

BOOST - we use Boost 1.55.0 but most recent versions are probably fine

####Installation
```
git clone https://github.com/sequencing/NxTrim.git
cd NxTrim
make
./nxtrim
```

You will also need to point the BOOST_ROOT environment variable at your boost installation eg.

``
export BOOST_ROOT=/your/boost/installation
``

if boost is installed globally then

``
export BOOST_ROOT=/usr/
``

should work

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

Trimming and assembly with [SPAdes](http://bioinf.spbau.ru/spades):
```
nxtrim -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -O sample --justmp
cat sample.mp.fastq.gz sample.unknown.fastq.gz > sample.allmp.fastq.gz
spades.py -k 21,33,55,77 -t 4 --hqmp1-12 sample.allmp.fastq.gz -o output_dir
```
Note we concatenate the unknown/mp libraries for SPAdes.  SPAdes versions>3.1.0 seem to dislike our virtual single/pe libraries hence we trim with the `--justmp` flag. This command is suitable for 2x151bp data, if you have 2x251bp then use `-k 21,33,55,77,127`.  

####Output:

The default behaviour expects raw fastq files from a Nextera Mate-Pair library kit in Reverse-Forward orientation.  Based on the location of the Nextera adapter sequence (if detected), nxtrim produces four different "virtual libraries":

* mp: read pairs that are large insert-size mate-pairs
* pe: read pairs that are short insert-sze paired-end reads
* se: single reads 
* unknown: a library of read-pairs that are mostly large-insert mate-pair, but possibly contain a small proportion of paired end contaminants

####Options:

The trimmer will reverse-complement the reads such that the resulting libraries will be in Forward-Reverse orientation, this reverse-complementing can be disabled via the --norc flag.

If you wish to generate pure mate-pair libraries (say for scaffolding), you can use the --justmp flag.  This will only generate the unknown and mp libraries.  Reads with an adapter occurring < minlength bp before the start will be completely N masked.

If you wish to preserve mate-pair libraries whenever possible, the --preservemp flag may be useful.  This will always keep the mate-pair library *unless* a read generated would be <minlength, in which case it will generate a PE.

You can trade specificity/sensitivity of adapter detection with the --similarity flag (1 - proportion of bp differences allowed for match) and the --minoverlap flag (minimum #bp considered on the ends of reads to match with the Nextera adapter).  The defaults were well suited to bacteria in our testing.

####Example data:

https://basespace.illumina.com/s/TXv32Ve6wTl9

Free registration required.

####References:

http://res.illumina.com/documents/products/technotes/technote_nextera_matepair_data_processing.pdf

http://res.illumina.com/documents/products/appnotes/appnote-nextera-mate-pair-bacteria.pdf
