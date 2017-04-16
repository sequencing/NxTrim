import sys

def N50(dna):
    lens = [len(val) for val in dna]
    lens.sort()
    thr = sum(lens)/2
    n50 = 0
    for i,l in enumerate(lens):
        n50+=l
        if n50>thr:
            return l

def read_assembly(fname,minlen):
    contigs = []
    scaffolds = []

    count=0
    s = []

    def parse(s):
        tmp = "".join(s)
#        print len(tmp)
        scaffolds.append(tmp)
        for c in [val for val in tmp.split("N") if len(val)>0]:
            contigs.append(c)
    infile = open(fname)
    line = infile.next()
    for linenum,line in enumerate(infile):
        if line[0]==">":
            if linenum>0:
                parse(s)
            s = []
            hdr = line
        else:
            count += len(line)
            s.append(line.strip())
    parse(s)

    return [val for val in contigs if len(val)>=minlen],[val for val in scaffolds if len(val)>=minlen]

if __name__ == "__main__":

    try: MINLEN = int(sys.argv[2])
    except: MINLEN=0

    contigs,scaffs = read_assembly(sys.argv[1],MINLEN)
    contig_lengths = [len(val) for val in contigs]
    scaff_lengths = [len(val) for val in scaffs]
    print "Assembly length\t%d"%sum(contig_lengths)
    print "#contigs\t%d"%len(contigs)
    print "Contig N50\t%d"%N50(contigs)
