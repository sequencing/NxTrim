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
git clone blah
cd nxtrim
make 

##Brief introduction to Nextera Mate Pair Libraries

The Nextera mate pair library involves circularising large (up to 12kb with a median of ~3kb) DNA fragments with the junction adapter:

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
##Adapter trimming logic

Clearly the artificial DNA sequence in the adapter needs to be removed before any analysis.  We now describe our adapter trimming routine, which attempts to maximise the amount of genomic DNA retained.  We achieve this by creating four "virtual" libraries:

	1. MP: Known mate-pair (long insert, RF orientation)
	2. UNKNOWN: Unknown (but likely mate-pair with long insert RF orientation)
	3. PE: Paired-end (short insert, FR orientation)
	4. SE: Single read (read has no mate)

Let (a1,b1) be the positions of the first and last base of the adapter (if found) in read one, and (a2,b2) be these positions for read two and L is the read length (same for both reads).  We set (a1,a2)=(L,L) and (a1,a2)=(L,L) when an adapter is not detected in respective reads.  We used R1[i,j] and R2[i,j] to refer to the substring from position i to j in reads 1 and two respectively.  We also set the variable MINLEN (default=25), which is the minimum read length (post-trimming) to report.

The virtual library construction logic is then:

```
    IF(a1>L AND a2>L) {
    	  Reverse-complement the pair and place in UNKNOWN
    }
    ELSE {
	IF(a1>MINLEN AND a1<L AND a2<MINLEN)
	     Place R1[0,a1] in SE. Discard R2.

	ELSE IF(a2>MINLEN AND a2<L AND a1<MINLEN)
	     Place R1[0,a1] in SE Discard R2.

	ELSE IF( (a1<L AND a2<L) OR ( a1<L AND b1>(L-MINLEN) ) OR ( a1<L AND b1>(L-MINLEN) ) ) 
	     Place reverse complement of the pair (R1[0,L],R2[0,L]) in MP

	ELSE IF b1<L AND a2==L {
	     IF a1<MINLEN
		Place (R1[b1,L],R2[0,L]) in PE
	     ELSE IF (L-b1)>a1 {
	     	Place (R1[b1,L],R2[0,L]) in PE
		Place R1[0,a1] in SE
	     }
	     ELSE {
		IF (L-b1)>MINLEN
		   Place R1[b1,L] in SE

	     	Place reverse complement of the pair (R1[0,b1],R2[0,L]) in MP
             }
	}		   	
	ELSE IF (b2<L AND a1==L) {
	     IF a2<MINLEN
		Place (R1[0,L],R2[b2,L]) in PE
	     ELSE IF (L-b2)>a2 {
	     	Place (R1[0,L],R2[b2,L]) in PE
		Place R2[0,a2] in SE
	     }
	     ELSE {
		IF (L-b2)>MINLEN
		   Place R2[b2,L] in SE

	     	Place reverse complement of the pair (R1[0,L],R2[0,b2]) in MP
             }	     
	}
```

##Adapter detection rules
We now describe the routine used to detect the presence/location of the adapter sequence in each read. There are two user defined variables involved, the similarity measure SIM (default=0.85) and the minimum string length to consider MINOVERLAP (default=12).  The similarity measure allows for adapter sequence with some amount of error to be matched, the MINOVERLAP variable allows us to consider partial adapter sequence on the very  beginning or end of a read.  The user also defines the distance(s1,s2) function as the Hamming (default) or Levenshtein distance between strings s1 and s2.  We describe the process for the DNA sequence R which could be either of the members of a read pair.

```

function detectAdapter(R,ADAPTER,SIM,distance) {
    start = MINOVERLAP - len(ADAPTER)
    stop = L - MINOVERLAP
    min_distance = len(ADAPTER)
    min_index = L
    distance_threshold = ceiling( (1 - SIM) * comparison_length)

    FOR i IN [start...stop] {
    	read_start = max(0,i)
	read_end = min(L,i)
	IF(i<0) {
		adapter_start = -i
        } ELSE {
                adapter_start = 0
        }
	IF( i>(L-len(ADAPTER) ) {
		adapter_end =  (L-i)
        } ELSE {
		adapter_end = len(adapter1)
        }		
    	comparison_length = adapter_end - adapter_start	
	d = distance(R[read_start:read_end],ADAPTER[adapter_start:adapter_end])
	IF( d<min_distance AND d<distance_threshold) {
	   min_distance = d
	   min_index = i
        }

   RETURN(min_index)
}

ADAPTER1 = CTGTCTCTTATACACATCT
ADAPTER2 = AGATGTGTATAAGAGACAG

a1 = detectAdapter(R1,ADAPTER1,SIM,distance)
IF(a1<L) {
   a1 = detectAdapter(R1,ADAPTER2,SIM,distance)
   IF(a1<L) {
      a1-=len(ADAPTER1)
   }
}
b1=a1+len(ADAPTER1)+len(ADAPTER1)

a2 = detectAdapter(R1,ADAPTER1,SIM,distance)
IF(a2<L) {
   a2 = detectAdapter(R1,ADAPTER2,SIM,distance)
   IF(a2<L) {
      a2-=len(ADAPTER1)
   }
}
b2=a2+len(ADAPTER1)+len(ADAPTER1)

```
We search for both the adapter sequences separately for two reasons;

   1. Occasionally only one of the adapters is present
   2. One of the adapters may have indel errors, but both having indel errors is unlikely. Hence we can use cheaper Hamming distance calculations for comparisons if we check them separately

##References:

http://res.illumina.com/documents/products/technotes/technote_nextera_matepair_data_processing.pdf

http://res.illumina.com/documents/products/appnotes/appnote-nextera-mate-pair-bacteria.pdf
