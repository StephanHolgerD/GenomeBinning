import pysam
#import pandas as pd
import sys
from multiprocessing import Pool
import argparse
import numpy as np
from tqdm import tqdm
from time import time
#!/usr/bin/env python3

t=time()
parser = argparse.ArgumentParser(description='Binning of mappings')
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-b','--bam', help='Pos. sorted and indexed bam file', required=True)
requiredNamed.add_argument('-o','--out', help='out file', required=True)

optArguments = parser.add_argument_group('optional arguments')
optArguments.add_argument('-c','--chromosome', help='name of chromosome/contig', required=False,default='all')
optArguments.add_argument('--bed', help='three column bed file for ROI', required=False,default=False)

optArguments.add_argument('-s','--binsize', help='size of bins', required=False,default=5000, type=int)
optArguments.add_argument('-p','--processes',default=1, help="number of cpu's  to run in paralell, ROI <1000 will always use 1 core",type=int)
#optArguments.add_argument('--tagKey',default=None, help="key for tag of interest",type=str)
#optArguments.add_argument('--tagString',default=None, help="string contained in tag of interest",type=str)







class BinChrom():
    def __init__(self,bam,binsize=100_000,chrom='all',tagKey=None,tagOfInterest=None,processes=12,bed=False):
        self.chrom=chrom
        self.bam=bam
        self.binsize=binsize
        self.bed=bed
        self.GetChromLenght()
        #self.tag=False
        #if tagKey is not None:
    #        self.tag=True
    #    self.tagKey=tagKey
    #    self.tagOfInterest=tagOfInterest
        self.processes=processes
    def GetChromLenght(self):
        self.ChrLen={}
        with pysam.AlignmentFile(self.bam) as f:
            if self.chrom != 'all':
                self.ChrLen[self.chrom]=f.get_reference_length(self.chrom)
            else:
                i=f.get_index_statistics()
                for c in i:
                    self.ChrLen[c.contig]=f.get_reference_length(c.contig)

    def chunks(self, chrom):
        chunkss=[]
        n=self.binsize

        if self.bed:
            with open(self.bed) as f:
                ROI=[[x.split('\t')[0],int(x.split('\t')[1])-1,int(x.split('\t')[2].rstrip())] for x in f if x.split('\t')[0] == chrom]  


            for roi in ROI:
                start=roi[1]
                end=roi[2]
                lst=range(start,end)
                for nn,i in enumerate(range(start,start+len(lst), n)):

                    i=i-start
                    if nn==0:
                        chunkss.append([roi[0],lst[i:i + n][0],lst[i:i + n][-1]])
                    else:
                        chunkss.append([roi[0],lst[i-1:i + n][0],lst[i:i + n][-1]])

            return chunkss

        else:
            lst=range(self.ChrLen[chrom])
            for i in range(0, len(lst), n):
                chunkss.append([chrom,lst[i:i + n][0],lst[i:i + n][-1]])
            return chunkss

    def FetchReads(self,chrom,start,end):
        cov=0
        r=0
        s=[]
        e=[]
        with pysam.AlignmentFile(self.bam) as f:
            for read in f.fetch(chrom,start,end):
                if not read.is_unmapped:
                    #print(read)
                    s.append(read.reference_start)
                    e.append(read.reference_end)
                    r=r+1
        s=[x if x >= start else start for x in s ]
        e=[x if x <= end else end for x in e ]
        s=np.array(s)
        e=np.array(e)
        nts=sum(e-s)
        cov=nts/(end-start)
        return [chrom,start,end,r,cov]

    #    else:
    #        cov=0
    #        r=0
    #        with pysam.AlignmentFile(self.bam) as f:
    #            for read in f.fetch(chrom,start,end):
    #                if not read.is_unmapped:
    #                    x=[x for x in read.tags if self.tagKey in x]
    #                    if len(x)<1:
    #                        continue
    #                    x=x[0][-1]
#
#                        if x==self.tagOfInterest:
#                            r=r+1
#                            cov=cov+read.rlen
#            cov=cov/(end-start)
#            return [start,end,cov,reads]

    def StartCalc(self):
        o = open(args.out,'w')
        o.write('contig,start,end,reads,coverage\n')

        for c in self.ChrLen:
            print(c)
            chunk=self.chunks(c)
            results=[]
            with Pool(processes=self.processes) as pool:
                jobs=[pool.apply_async(self.FetchReads,args=(i[0],i[1],i[2])) for i in chunk]

                for job in tqdm(jobs):
                    results.append(job.get())
            for r in results:
                o.write(f'{str(r[0])},{str(r[1])},{str(r[2])},{str(r[3])},{str(r[4])}\n')


if __name__ == '__main__':
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args=parser.parse_args()
    #cc=BinChrom(args.bam,tagKey=args.tagKey,tagOfInterest=args.tagString,binsize=args.binsize,chrom=args.chromosome,processes=args.threads)
    cc=BinChrom(args.bam,binsize=args.binsize,chrom=args.chromosome,processes=args.processes, bed=args.bed)

    cc.StartCalc()
   
print(time()-t)
