nxtrim: Software to remove Nextera Mate Pair adapters and categorise reads according to the orientation implied by the adapter location.

Copyright (c) 2014, Illumina, Inc.
All rights reserved.

See LICENSE.txt for details (BSD License).

####Dependencies

BOOST


####Installation
```
git clone git@github.com:sequencing/NxTrim.git
cd NxTrim
make
```
####Usage
Trim the data:
```
nxtrim -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -O sample --rc
```

Assemble with Velvet:
```
velveth output_dir 55 -short -fastq.gz sample.se.fastq.gz -shortPaired2 -fastq.gz sample.pe.fastq.gz -shortPaired3 -fastq.gz sample.mp.fastq.gz -shortPaired4 -fastq.gz sample.unknown.fastq.gz
velvetg output_dir -exp_cov auto -cov_cutoff auto -shortMatePaired4 yes
```

Assemble with SPAdes:
```
cat sample.mp.fastq.gz sample.unknown.fastq.gz > sample.allmp.fastq.gz
spades.py -k 21,33,55,77 -t 4 --careful --pe1-s sample.se.fastq.gz --pe2-12 sample.pe.fastq.gz --hqmp3-12 sample.allmp.fastq.gz --hqmp3-fr --careful -o output_dir
```
Note we concatenate the unknown/mp libraries for SPAdes.  This command is suitable for 2x151bp data, if you have 2x251bp then use `-k 21,33,55,77,127`

####References:

http://res.illumina.com/documents/products/technotes/technote_nextera_matepair_data_processing.pdf

http://res.illumina.com/documents/products/appnotes/appnote-nextera-mate-pair-bacteria.pdf

####Example data:

https://basespace.illumina.com/s/TXv32Ve6wTl9

Free registration required.
