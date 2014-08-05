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
./nxtrim -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -O sample --rc
```
Assemble with velvet:
```
velveth output_dir 55 -short -fastq.gz sample.se.fastq.gz -shortPaired2 -fastq.gz sample.pe.fastq.gz -shortPaired3 -fastq.gz sample.mp.fastq.gz -shortPaired4 -fastq.gz sample.unknown.fastq.gz

velvetg output_dir -exp_cov auto -cov_cutoff auto -shortMatePaired4 yes
```

####References:

http://res.illumina.com/documents/products/technotes/technote_nextera_matepair_data_processing.pdf

http://res.illumina.com/documents/products/appnotes/appnote-nextera-mate-pair-bacteria.pdf

####Example data:

https://basespace.illumina.com/s/TXv32Ve6wTl9

Free registration required.
