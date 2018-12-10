from dataclasses import dataclass, field
import collections
import numpy as np
from typing import List, Dict
from argparse import ArgumentParser, RawTextHelpFormatter
import re
import gzip
from Bio.Seq import Seq
import pandas as pd
import copy

@dataclass
class BWT():
    T: str
    inds: List[str] = field(default_factory=list)
    charfreq: Dict[str, int] = field(default=None)
    charspace: List[str] = field(default_factory=list)
    start: int = 29130001 # include start coordinate, 1-base
    stop: int = 29140000  # include stop coordinate, 1-base
    offset: int = 29130000
    
    
    def vcl(self):
        """ Get alternative variant call from vcf files.
            Store the variants into snps dictionary."""
        self.snps = dict()
        with gzip.open('../raw_data/ALL.chr22.SNPs.10k.vcf.gz', 'rb') as handle:
            dfSNP = pd.read_csv(handle, sep='\t')
            dfTmp = dfSNP[(dfSNP['POS']>self.start) & (dfSNP['POS']<self.stop)]
            dfTmp = dfTmp[dfTmp['ALT'].apply(lambda x: len(x)==1)].reset_index(drop=True)
            for dfi, row in dfTmp.iterrows():
                self.snps[row['POS']] = row['ALT'][0]
            
        print("SNP list complete!")
    
    def genfirst(self):
        """ Generate first column called F array in LF(Last First) algorithm.
            The F array is simplified into a collection of sorted alphabet with frequency occurred in T."""
        self.charfreq = collections.Counter(self.T)
        self.charspace = sorted(self.charfreq.keys())
    
    def encodeT(self, i):
        refchar = self.T[i-1]
        nuc_ar = np.array(['A', 'C', 'G', 'T'])
        if i in self.snps.keys():
            altchar = self.snps[i]
            return (nuc_ar == refchar) | (nuc_ar == altchar)
        else:
            return nuc_ar == refchar
    
    def gencp(self):
        """ Generate checkpoints Tally array for traversal"""
        chars = np.array(self.charspace[1:], dtype=np.string_)
        freqs = np.zeros(chars.size)
        self.cp = []
        for c in self.bwstr:
            freqs = freqs + (chars == c.encode())
            self.cp.append(freqs)
        self.cp = np.array(self.cp)
        
    def toSA(self):
        """ Create suffix array and return the index sorted by lexicographical order"""
        lT = len(self.T)
        sat = sorted([(self.T[i:lT], i) for i in range(lT)])
        self.sufidx = list(map(lambda x: x[1], sat))
        return self.sufidx
    
    def satobwt(self):
        """ Return BWT string using SA"""
        bw = []
        self.bwstr = ''
        for i in self.toSA():
            if i==0:
                bw.append(np.array([0,0,0,0], dtype=np.bool))
                self.bwstr+='$'
            else:
                bw.append(self.encodeT(i))
                self.bwstr+= self.T[i-1]
        return np.array(bw)
    
    def query(self, q):
        self.inds = self.inds.replace('\'', '').replace('[', '').replace(']', '').replace(' ', '').split(',')
        
        q = q.replace(' ', '').upper() 
        rev_q = str(Seq(q).reverse_complement())
        re_f = self._query(q, '+')
        re_r = self._query(rev_q, '-')
        re = re_f + re_r
        if re:
            print("Pattern found:\n%s" %re)
        else:
            print('No pattern matched.')
            
    def _query(self, q, strand):
        result = []
        chars = np.array(self.charspace[1:], dtype=np.string_)
        self.sufidx = list(self.sufidx)
        frestidx, scores, curstr = self._finitidx(q[-1])
        #print("Initial first index: %s" % str(list(frestidx)))
        if frestidx: # if there's any initial index left
            for c in q[:-1][::-1]:
                #print('current checking c %s' % c)
                tmpidx = []
                for i in frestidx:
                    
                    if self.bwstr[i]==c:
                        #print("bwt %d is %s" % (i, self.bwstr[i]))
                        crank = int(self.cp[i, np.argwhere(chars==c.encode()).flatten()[0]])
                        #print("crank is %d" % crank)
                        nextidx = self._findidx(c, crank)
                        #print('nextidx = %d' % nextidx)
                        tmpidx.append(nextidx)
                        
                frestidx = tmpidx
                #print(frestidx)
            foundidx = self._resolve(frestidx)
            if foundidx:
                gRNA = self.T[foundidx[0]:(foundidx[0]+ len(q))]
                if strand =='-':
                    gRNA = str(Seq(gRNA).reverse_complement())
                result = []
                for idx in foundidx:
                    idx = (idx+self.offset)
                    ar, dfSNP = self._containsnp(idx, q)
                    table = ['gRNA', 'start', 'stop', 'SNP ID', 'strand', 'MIT score', 'groups']
                    groups = dict()
                    groups[q] = self.inds
                    if ar:
                        
                        for dfid, row in dfSNP.iterrows():
                            newg = dict()
                            newg = copy.deepcopy(groups)
                            for k in groups.keys():
                                
                                #print("k is %s"% groups[k])
                                dfTmp = dfSNP[groups[k]]
                                ref_grna = k
                                alt_groups = dfTmp.columns[row[newg[k]] != '0|0']
                                ref_groups = dfTmp.columns[row[newg[k]] == '0|0']
                                if alt_groups.any():
                                    #print("Who is in the group for alternative: %s" % alt_groups)
                                    grna_pos = row['POS'] - idx -1
                                    alt_grna = q[:grna_pos] + row['ALT'] + q[(grna_pos+1):]
                                    #print("alt_grna is %s; ALT is %s; REF is %s; POS is %s" % (alt_grna, row['ALT'], row['REF'], row['POS']))
                                    if not alt_grna in newg:
                                        newg.update({alt_grna:list(alt_groups)}) 
                                newg.update({ref_grna:list(ref_groups)})
                                
                            
                            groups = copy.deepcopy(newg)
                        
                        for k in groups:
                            result.append({'gRNA': gRNA,
                                           'Target': k,
                                           'Cleavage efficiency': self._score(q,k, strand),
                                           'Start': (idx+1),
                                           'Stop': (idx+len(q)),
                                           'Strand': strand,
                                           'Sample ID': groups[k],
                                           })
                            
                            
                    else:
                        print('No SNP found in gRNA targeting region.')
                return result
            else:
                return result
        else:
            return result
    
    def _score(self, q, k, strand):
        self.penalties = np.array([0, 0, 0.014, 0, 0, 0.395, 0.317, 0,
                                   0.389, 0.079, 0.445, 0.508, 0.613,
                                   0.851, 0.732, 0.828, 0.615, 0.804,
                                   0.685, 0.583, 0, 0.75, 1])
        if strand == '-':
            self.penalties = self.penalties[::-1]
        score = float(1)
        try:
            for qc, kc, pc in zip(q, k, self.penalties[-len(q):]):
                if qc!=kc:
                    score = score * (1-pc)
            return score
        except IndexError:
            print("Length of gRNA is longer than the size of penalty matrix.")
                    
    
    def _containsnp(self, idx, q):
        try:
            
            with gzip.open('../raw_data/ALL.chr22.SNPs.10k.vcf.gz', 'rb') as handle:
                dfSNP = pd.read_csv(handle, sep='\t')
                dfTmp = dfSNP[(dfSNP['POS']>idx) & (dfSNP['POS']<(idx+len(q)))]
                dfTmp = dfTmp[dfTmp['ALT'].str.len() ==1].reset_index(drop=True)
            print("Matched position:%d-%d" % (idx, (idx+len(q))))

            return True, dfTmp
        
        except:
            return False, _
    
    def _findidx(self, c, crank):
        idx=0
        for s in self.charspace:
            if s == c:
                idx = idx + crank - 1
                return idx
            idx += self.charfreq[s]
        return idx 
        
    def _finitidx(self, c):
        idx=0; end=0
        for s in self.charspace:
            if s == c:
                end = idx+self.charfreq[s]
                break
            else:
                idx += self.charfreq[s]
        fidxs = []
        scores = []
        curstr = []
        for i in range(idx, end):
            fidxs.append(i)
            scores.append([1])
            curstr.append([c])
        return fidxs, scores, curstr
    
    def _resolve(self, lidx):
        #print("Resolving %s" % (lidx))
        found = []
        chars = np.array(self.charspace[1:], dtype=np.string_)
        for i in lidx:
            found.append(int(self.sufidx[i]))
            #print("CURRENT FOUND %s" % found)
        return found
    
    def main(self):
        self.T = self.T.upper() + '$'
        self.vcl()
        self.genfirst()
        self.bw = self.satobwt()
        self.gencp()
        print("BWT complete!")
        #print(self.bw)

def main(args):
    att = ['txt', 'fa', 'fasta', 'gz']
    if args.i.split('.')[-1] in att:
        if args.i.split('.')[-1] == 'gz':
            t = ''
            with gzip.open(args.i, 'rb') as handle:
                fl = handle.readline().decode().strip('\w').rstrip('\n')
                if re.search('^>.*', fl):
                    rl = handle.readlines()
                    t = ''.join([l.decode().strip('\w').rstrip('\n') for l in rl])

                elif re.search('^[ATCGNatcgn]+$', fl):
                    rl = handle.readlines()
                    t = fl + ''.join([l.decode().strip('\w').rstrip('\n') for l in rl])
                else:
                    print("File format has something wrong.")
                    exit(1)
            print("Reference sequence file loading %s" % args.i)
    elif re.search('^[ATCGNatcgn]+$', args.i):
        t = args.i.replace(' ', '').rstrip('\n')
        print("Reference sequence string input: %s" % args.i)
        
    bb = BWT(T=t, inds=args.g)
    bb.main()
    bb.query(args.q)
    
        
if __name__ == '__main__':
    
    parser = ArgumentParser(prog='bwt', 
                            description='BWT for CRISPR/Cas9',
                            formatter_class=RawTextHelpFormatter)
    
    parser.add_argument("-i", "--file", help="[FILE] or str\n[FILE] accepts attributes ['txt', 'fa', 'fasta']",
                        required=True, dest='i', metavar='')
    parser.add_argument("-q", "--query", help="Sequence pattern for search in given reference genome file.",
                        required=True, dest='q', metavar='')
    parser.add_argument("-g", "--group", help="List of individuals selected",
                        required=True, dest='g', metavar='')
    #try:
    args = parser.parse_args()
    main(args)
        
    #except:
    #    print("Please enter a valid string or a file with DNA sequence.")
    #    parser.print_help()
    #    exit(0)