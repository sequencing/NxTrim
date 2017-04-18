### Adapter matching:

We implement a very simple approach. Each read is searched for the junction adapter (`CTGTCTCTTATACACATCT`) and its reverse complement (`AGATGTGTATAAGAGACAG`) based on Hamming distance, the entire junction being `CTGTCTCTTATACACATCT`+`AGATGTGTATAAGAGACAG`. If a match to either side of the junction is found then the entire 38bp juncton is clipped. This means the algorithm can tolerate an indel error in one of the sides of the junction, but not both. Optionally Smith-Waterman alignment can be turned on with `-w` which resolves this issue, although the code is less tested.

eg. this will be caught:

```
CTGTCTCT-TA-TACACATCTAGATGTGTATAAGAGACAG
```

but not this:

```
CTGTCTCT-TA-TACACATCTAGATGT-GTATAAGAGACAG
```
