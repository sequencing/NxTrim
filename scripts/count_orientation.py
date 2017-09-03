import sys,pysam,itertools,time,collections


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='count orientation of read pairs (for unsorted bam)')
    parser.add_argument('bam', metavar='bam', type=str, help='bam')
    args = parser.parse_args()
    sam = pysam.Samfile(args.bam, "rb" )
    out_fr = pysam.Samfile(args.bam[:-3]+"fr.bam", "wb", header = sam.header)
    out_rf = pysam.Samfile(args.bam[:-3]+"rf.bam", "wb", header = sam.header)
    L = 151
    d = collections.defaultdict(int)
    thr=20000 #how much to clip off contig ends?
    
    for read in sam:
        not_on_end = read.pos>thr and read.pos < (sam.lengths[read.tid] - thr) and read.pnext>thr and read.pnext< (sam.lengths[read.tid] - thr) 
        not_too_small = (max(read.pos+L,read.pnext+L) - min(read.pos,read.pnext))> (L)
        not_identical = read.pos!=read.pnext

        if not read.is_unmapped and not read.mate_is_unmapped and not_identical:            
            if not_on_end and not_too_small:
                if read.pos < read.pnext:
                    if (read.flag & 0x10):
                        k1 = 'R'
                    else: 
                        k1 = 'F'

                    if (read.flag & 0x20):
                        k2='R'
                    else:  
                        k2='F'
                    o=k1+k2
                    d[k1+k2] += 1
                else:
                    if (read.flag & 0x10):
                        k1 = 'R'
                    else: 
                        k1 = 'F'

                    if (read.flag & 0x20):
                        k2='R'
                    else:  
                        k2='F'
                    o=k2+k1

                if o=='RF':
                    out_rf.write((read))
                if o=='FR':
                    out_fr.write((read))
                

    for k in ['FR','RF','FF','RR']:
        print k,d[k]
