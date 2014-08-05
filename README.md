nxtrim: Software to remove Nextera Mate Pair adapters and categorise reads according to the orientation implied by the adapter location.

##LICENSE

Copyright (c) 2014, Illumina, Inc.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

##DEPENDENCES

BOOST


##INSTALLATION
```
git clone git@github.com:sequencing/NxTrim.git
cd NxTrim
make
```
##USAGE

```
./nxtrim -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -O sample --rc
```

##BACKGROUND

The Nextera mate pair libraries involve circularising large (up to 12kb with a median of ~4kb) DNA fragments with the junction adapter:

```
    CTGTCTCTTATACACATCT+AGATGTGTATAAGAGACAG
```

this circularised DNA is then fragmented again into smaller typical PE sized chunks (several hundred bp).  The chunks containing this adapter are then sequenced.  DNA on either side of the adapter sequence will be in Reverse-Forward orientation with a mate-pair insert size.  The can adapter occur (roughly) uniformly across the shorter fragment. We now enumerate the possible outcomes of the adapter location, note we use arbitrarily short read lengths for readability, in practice 2 x 151bp (or larger) lengths can be expected.


For example we will commonly observe read-pairs containing the partial or whole adapter sequence (we use X to denote true DNA from the left of the adapter and Y for the right):

```
R1--------------------------------------------------->
  XXXXXXXXXXXXXXXXXXXXXXXXCTGTCTCTTATACACATCTAGATGTGTATAAGAGACAGYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY
                                                                                                  <---------------------------------------------------R2
```

There are cases where the adapter is not observed either because

     1. the reads were not long enough to reach the adapter:
```
R1--------------------------------------------------->
  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXCTGTCTCTTATACACATCTAGATGTGTATAAGAGACAGYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY
                                                                                                  <---------------------------------------------------R2
```
     2. the read is a PE contaminant and had no adapter:
```
R1--------------------------------------------------->
  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                                                                                                  <---------------------------------------------------R2
```
the first case is far more frequent, but there can be a small amount of contamination from PE reads.

Another thing to consider is that the adapter may occur very early in one of the reads, meaning the preceding DNA (if any) is far to short to be of use for assembly:
```
R1--------------------------------------------------->
  XXCTGTCTCTTATACACATCTAGATGTGTATAAGAGACAGYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY
                                                                                                  <---------------------------------------------------R2
```

Clearly the artificial DNA sequence in the adapter needs to be removed before any analysis.  We now describe our adapter trimming routine, which attempts to maximise the amount of genomic DNA retained.  We achieve this by creating four "virtual" libraries:

	1. MP: Known mate-pair (long insert, RF orientation)
	2. UNKNOWN: Unknown (but likely mate-pair with long insert RF orientation)
	3. PE: Paired-end (short insert, FR orientation)
	4. SE: Single read (read has no mate)

##References:

http://res.illumina.com/documents/products/technotes/technote_nextera_matepair_data_processing.pdf

http://res.illumina.com/documents/products/appnotes/appnote-nextera-mate-pair-bacteria.pdf
