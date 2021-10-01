import pysam
#import pandas as pd
import sys
from multiprocessing import Pool
import argparse



#!/usr/bin/env python3

parser = argparse.ArgumentParser(description='Binning of mappings')
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-b','--bam', help='Pos. sorted and indexed bam file', required=True)
requiredNamed.add_argument('-o','--out', help='out file', required=True)

optArguments = parser.add_argument_group('optional arguments')
optArguments.add_argument('-c','--chromosome', help='name of chromosome/contig', required=False,default='1')
optArguments.add_argument('-s','--binsize', help='size of bins', required=False,default=5000)
optArguments.add_argument('--threads',default=1, help="number of cpu's  to run in paralell, ROI <1000 will always use 1 core",type=int)
optArguments.add_argument('--tagKey',default=None, help="key for tag of interest",type=str)
optArguments.add_argument('--tagString',default=None, help="string contained in tag of interest",type=str)







class BinChrom():
    def __init__(self,bam,binsize=100_000,chrom='1',tagKey=None,tagOfInterest=None,processes=12):
        self.chrom=chrom
        self.bam=bam
        self.binsize=binsize
        self.GetChromLenght()
        self.tag=False
        if tagKey is not None:
            self.tag=True
        self.tagKey=tagKey
        self.tagOfInterest=tagOfInterest
        self.processes=processes
    def GetChromLenght(self):
        with pysam.AlignmentFile(self.bam) as f:
            self.ChrLen=f.get_reference_length(self.chrom)
            i=f.get_index_statistics()
            i=[x for x in i if x.contig == self.chrom]
            ii=i[0].mapped
            ii=ii+i[0].unmapped
            self.Reads=ii

    def chunks(self):
        chunkss=[]
        lst=range(self.ChrLen)
        n=self.binsize
        for i in range(0, len(lst), n):
            chunkss.append([lst[i:i + n][0],lst[i:i + n][-1]])
        return chunkss

    def FetchReads(self,start,end):
        if not self.tag:
            cov=0
            with pysam.AlignmentFile(self.bam) as f:
                for read in f.fetch(self.chrom,start,end):
                    if not read.is_unmapped:
                        cov=cov+read.rlen
            cov=cov/(end-start)
            return [start,end,cov]

        else:
            cov=0
            with pysam.AlignmentFile(self.bam) as f:
                for read in f.fetch(self.chrom,start,end):
                    if not read.is_unmapped:
                        x=[x for x in read.tags if self.tagKey in x]
                        if len(x)<1:
                            continue
                        x=x[0][-1]

                        if x==self.tagOfInterest:
                            cov=cov+read.rlen
            cov=cov/(end-start)
            return [start,end,cov]

    def StartCalc(self):
        chunk=self.chunks()
        print('calc chunk')
        with Pool(processes=self.processes) as pool:
            results=pool.starmap(self.FetchReads,chunk)
        return results


if __name__ == '__main__':
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args=parser.parse_args()
    cc=BinChrom(args.bam,tagKey=args.tagKey,tagOfInterest=args.tagString,binsize=args.binsize,chrom=args.chromosome,processes=args.threads)
    print('start multi')
    results=cc.StartCalc()
    with open(args.out,'w') as o:
        o.write('s,e,c\n')
        for r in results:
            o.write(f'{str(r[0])},{str(r[1])},{str(r[2])}\n')
