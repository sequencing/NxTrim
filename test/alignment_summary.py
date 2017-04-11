import sys,pysam,collections


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='routine to ntrim adapter sequences from Nextera mate pair data')
    parser.add_argument('bam', metavar='bam', type=str, help='bam')
    parser.add_argument('-output', metavar='output', type=str, help='output prefix',default=None)
    args = parser.parse_args()
    sam = pysam.Samfile(args.bam, "rb" )

    if(args.output!=None):
        out_fr = pysam.Samfile(args.bam[:-3]+"fr.bam", "wb", header = sam.header)
        out_rf = pysam.Samfile(args.bam[:-3]+"rf.bam", "wb", header = sam.header)

    L = 151
    d = collections.defaultdict(int)
    isizes = collections.defaultdict(list)
    thr=20000 #how much to clip off contig ends?
    
    for read in sam:

        not_on_end = read.pos>thr and read.pos < (sam.lengths[read.tid] - thr) and read.pnext>thr and read.pnext< (sam.lengths[read.tid] - thr) 
        insert_size = max(read.pos+L,read.pnext+L) - min(read.pos,read.pnext)
        not_too_small = insert_size > L
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
                    isizes[k1+k2].append(insert_size)
#                    sys.stderr.write("%s %d\n"%(k1+k2,insert_size))
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
                if(args.output!=None):
                    if o=='RF':
                        out_rf.write((read))
                    if o=='FR':
                        out_fr.write((read))
                
    print "Orientation","Frequency".rjust(10," "),"Median insert size".rjust(20," ")
    for k in ['FR','RF','FF','RR']:
        x=isizes[k]
        x.sort()
        print k.ljust(9," "),("%d"%d[k]).rjust(12," "),("%d"%x[len(x)/2]).rjust(20," ")

    print "\nRF/(FR+RF) = %f"%(d['RF']/float(d['RF']+d['FR']))
