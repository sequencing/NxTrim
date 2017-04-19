## Adapter matching:

### Default behaviour:

We implement a very simple approach. Each read is searched for the junction adapter (`CTGTCTCTTATACACATCT`) and its reverse complement (`AGATGTGTATAAGAGACAG`) based on Hamming distance, the entire junction being `CTGTCTCTTATACACATCT`+`AGATGTGTATAAGAGACAG`. If a match to either side of the junction is found then the entire 38bp juncton is clipped. This means the algorithm can tolerate an indel error in one of the sides of the junction, but not both. 

eg. this will be caught:

```
CTGTCTCT-TA-TACACATCTAGATGTGTATAAGAGACAG
```

but not this:

```
CTGTCTCT-TA-TACACATCTAGATGT-GTATAAGAGACAG
```

### Aggressive mode

Turning on `--aggressive` will seek for adapters more aggressively. Rather than just check for the left/right side of the adapters with Hamming distance, it shreds the 38bp junction adapter into 19-mers and then checks if any of these are present in the reads. If any 19-mer match within the specified similarity is found, it will clip the read appropriately. This approach finds a very modest increase in adapters, but doesn't have a practical impact on assembly quality in my hands. Your mileage may vary.

Ideally we would use a seed-and-extend approach as implemented in tools such as [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic). The rate of missed adapters with `--aggressive` on was around 1 per 8000 read-pairs in one test on 2x151bp E. coli K-12 MG1655 run (missed adapters found with blast).

