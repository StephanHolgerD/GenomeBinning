import pysam
#import pandas as pd
import sys
from multiprocessing import Pool


class BinChrom():
    def __init__(self,bam,binsize=100_000,chrom='1',tag=False,tagKey='',tagOfInterest='',processes=12):
        self.chrom=chrom
        self.bam=bam
        self.binsize=binsize
        self.GetChromLenght()
        self.tag=tag
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
        with Pool(processes=self.processes) as pool:
            results=pool.starmap(self.FetchReads,chunk)
        return results


if __name__ == '__main__':
    cc=BinChrom(sys.argv[1],tag=True,tagKey='XS',tagOfInterest='Unassigned_NoFeatures',binsize=5_000)
    print('start multi')
    results=cc.StartCalc()
    with open(sys.argv[2],'w') as o:
        o.write('s,e,c\n')
        for r in results:
            o.write(f'{str(r[0])},{str(r[1])},{str(r[2])}\n')
