nxtrim: Software to remove Nextera Mate Pair adapters and categorise reads according to the orientation implied by the adapter location.  This software is not commercially supported.

Copyright (c) 2016, Illumina, Inc. All rights reserved.

This software is provided under the terms and conditions of the BSD 2-Clause License

You should have received a copy of the BSD 2-Clause License along with this program. If not, see https://github.com/sequencing/licenses/.

Some detailed assembly results for Nextera Mate-Pair data are available [here](https://github.com/sequencing/NxTrim/wiki/Bacterial-assembles-using-Nextera-Mate-pairs).

### Installation

```
git clone https://github.com/sequencing/NxTrim.git
cd NxTrim
make
./nxtrim
```

###Usage

Trimming and assembly with [Velvet](https://www.ebi.ac.uk/~zerbino/velvet/):

```
nxtrim -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -O sample 
velveth output_dir 55 -short -fastq.gz sample.se.fastq.gz -shortPaired2 -fastq.gz sample.pe.fastq.gz -shortPaired3 -fastq.gz sample.mp.fastq.gz -shortPaired4 -fastq.gz sample.unknown.fastq.gz
velvetg output_dir -exp_cov auto -cov_cutoff auto -shortMatePaired4 yes
```
the above approach corresponds to the results in the NxTrim publication.

Trimming and assembly with [SPAdes](http://bioinf.spbau.ru/spades):

```
nxtrim -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -O sample --justmp
cat sample.mp.fastq.gz sample.unknown.fastq.gz > sample.allmp.fastq.gz
spades.py -t 4 --hqmp1-12 sample.allmp.fastq.gz -o output_dir
```
We concatenate the unknown/mp libraries for SPAdes.  SPAdes versions>3.1.0 seems to perform better without our virtual single/pe libraries hence we trim with the `--justmp` flag. 

**Note:** We achieved good results using the above commands on the bacterial samples analysed in the NxTrim paper.  These had modest coverage (<50X).  If you have very high coverage samples, it might be preferable to not use the "unknown" library at all or just treat it as a single-ended library, this will remove the risk of PE contaminants causing problems.

Piping trimmed reads directly to an aligner:

```
nxtrim --stdout -1 EcMG1_ATGTCA_L001_R1_001.fastq.gz -2 EcMG1_ATGTCA_L001_R2_001.fastq.gz | bwa mem EcMG.fna -p - > out.sam
or
nxtrim --stdout-mp -1 EcMG1_ATGTCA_L001_R1_001.fastq.gz -2 EcMG1_ATGTCA_L001_R2_001.fastq.gz | bwa mem EcMG.fna -p - > out.sam
```
The first command pipes both unknown/MP reads to stdout, this is useful if you have a high quality reference to align to. The second only prints *known* MP reads, which is useful for scaffolding purposes.

###Output:

The default behaviour expects raw fastq files from a Nextera Mate-Pair library kit in Reverse-Forward orientation.  Based on the location of the Nextera adapter sequence (if detected), nxtrim produces four different "virtual libraries":

* mp: read pairs that are large insert-size mate-pairs
* pe: read pairs that are short insert-sze paired-end reads
* se: single reads 
* unknown: a library of read-pairs that are mostly large-insert mate-pair, but possibly contain a small proportion of paired end contaminants

###Options:

The trimmer will reverse-complement the reads such that the resulting libraries will be in Forward-Reverse (FR) orientation, if you wish to keep your reads as Reverse-Forward then use --rf flag.

If you wish to generate pure mate-pair libraries (say for scaffolding), you can use the --justmp flag.  This will only generate the unknown and mp libraries.  Reads with an adapter occurring < minlength bp before the start will be completely N masked.

If you wish to preserve mate-pair libraries whenever possible, the --preservemp flag may be useful.  This will always keep the mate-pair library *unless* a read generated would be <minlength, in which case it will generate a PE.

You can trade specificity/sensitivity of adapter detection with the --similarity flag (1 - proportion of bp differences allowed for match) and the --minoverlap flag (minimum #bp considered on the ends of reads to match with the Nextera adapter).  The defaults were well suited to bacteria in our testing.

###Example data:

https://basespace.illumina.com/s/TXv32Ve6wTl9

Free registration required.

###References:

http://res.illumina.com/documents/products/technotes/technote_nextera_matepair_data_processing.pdf

http://res.illumina.com/documents/products/appnotes/appnote-nextera-mate-pair-bacteria.pdf

[O’Connell, Jared, et al. "NxTrim: optimized trimming of Illumina mate pair reads." Bioinformatics 31.12 (2015): 2035-2037.](http://bioinformatics.oxfordjournals.org/content/31/12/2035.abstract)
