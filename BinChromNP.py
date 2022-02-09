import pysam
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

optArguments.add_argument('--chunksize', help='size of chunks', required=False,default=500_000, type=int)

optArguments.add_argument('-p','--processes',default=1, help="number of cpu's  to run in paralell, ROI <1000 will always use 1 core",type=int)


class BinChrom():
    def __init__(self,bam,binsize=100_000,chunksize=500_000,chrom='all',tagKey=None,tagOfInterest=None,processes=12,bed=False):
        self.chrom=chrom
        self.bam=bam
        self.binsize=binsize
        self.bed=bed
        self.chunksize=chunksize
        self.GetChromLenght()

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

    def chunks(self, chrom,size):
        chunkss=[]
        n=size

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


    def FetchReadsChunk(self,chrom,start,end):
        s=[]
        e=[]
        with pysam.AlignmentFile(self.bam) as f:
            for read in f.fetch(chrom,start,end):
                if not read.is_unmapped:
                    s.append(read.reference_start)
                    e.append(read.reference_end)
        s=np.array(s)
        e=np.array(e)

    def CalcBin(self,chrom,start,end,s_array,e_array):
        cov=0
        s1=np.where(((s_array <=start) & (e_array >= end)))
        cov=cov+self.binsize*len(s1[0])
        s2=np.where((s_array >=start)&(s_array <=end) & (e_array >= end))
        ss=s_array[s2[0]]
        cov=cov+(sum(end-ss))
        s3=np.where(((s_array <=start)&(e_array <=end) & (e_array >= start)))
        ee=e_array[s3[0]]
        cov=cov+(sum(ee-start))
        r=end-start
        if r != 0: 
            cov=cov/r
        else:
            cov=cov

        align=len(s1[0])+len(s2[0])+len(s3[0])
        return chrom,start,end,align,cov
    
    def CalcChunk(self,chrom,start,end):
        results=[]
        s=[]
        e=[]
        with pysam.AlignmentFile(self.bam) as f:
            for read in f.fetch('chr1',start, end):
                if not read.is_unmapped:
                    s.append(read.reference_start)
                    e.append(read.reference_end)
        s=np.array(s)
        e=np.array(e)
        lst=range(start,end+1)
        for a in range(0,len(lst),self.binsize):
            if len(e) == 0:
                chunks=lst[a:a+self.binsize][0]
                chunke=lst[a:a+self.binsize+1][-1]
                results.append([chrom,chunks,chunke,0,0])
            else:
                chunks=lst[a:a+self.binsize][0]
                chunke=lst[a:a+self.binsize+1][-1]
                contig,start,end,align,cov=self.CalcBin(chrom,chunks,chunke,s,e)
                results.append([contig,start,end,align,cov])
        return results

    def StartCalc(self):
        o = open(args.out,'w')
        o.write('contig,start,end,reads,coverage\n')
        for c in self.ChrLen:
            print(c)
            results=[]
            chunk=self.chunks(c, self.chunksize)
            with Pool(processes=self.processes) as pool:
                jobs=[pool.apply_async(self.CalcChunk,args=(i[0],i[1],i[2])) for i in chunk]
                for job in tqdm(jobs):
                    results.append(job.get())
            for _ in results:
                for r in _:
                    o.write(f'{str(r[0])},{str(r[1])},{str(r[2])},{str(r[3])},{str(r[4])}\n')
        o.close()    



if __name__ == '__main__':
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args=parser.parse_args()
    cc=BinChrom(args.bam,binsize=args.binsize,chrom=args.chromosome,processes=args.processes, bed=args.bed, chunksize=args.chunksize)

    cc.StartCalc()
   
print(time()-t)
