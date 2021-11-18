import glob
import os
import argparse 
import errno

import subprocess
from subprocess import call
from multiprocessing import Pool

import pandas as pd
import numpy


parser = argparse.ArgumentParser()
parser.add_argument('-m', '--raw_maf_dir')
parser.add_argument('-t', '--mafIndex_dir', required=True)
parser.add_argument('-ch', '--chromSizes', required=True) 
parser.add_argument('-n', '--dirname', required=True)
parser.add_argument('-g', '--global_coo_file', required=True)

args = parser.parse_args()

# global arguments - paths to files/tools
raw_maf_dir = args.raw_maf_dir
mafIndex_dir = args.mafIndex_dir
chromSizes = args.chromSizes
dirname = args.dirname
global_coo_file = args.global_coo_file
outdir_path = os.path.join(raw_maf_dir, dirname) 

print ('path to dir where to,  download maf, unpack them, get .bb index file: ', raw_maf_dir)
print ('path to mafIndex and mafExtract tools: ', mafIndex_dir)
print ('path to chromSizes file: ', chromSizes)
print ('path to older dir including bed, maf, fasta result files: ', dirname)
print ('path to global coo file: ', global_coo_file)

# create dir structure 
print ('__________________________')
print ('creating dir structure')

try:
    os.makedirs(os.path.join(raw_maf_dir, dirname))     
except Exception as e:
    print (e)
    
    
    
try:
    for chrom in ['chr'+str(i) for i in range(1, 23)]+['chrY', 'chrX']:
        print (os.path.join(outdir_path, chrom))
        os.makedirs(os.path.join(outdir_path, chrom))      
except Exception as e:
    print (e)        
    
    
    
try:
    for chrom in ['chr'+str(i) for i in range(1, 23)]+['chrY', 'chrX']:
        for extension in ['bed', 'maf', 'maf_m', 'fasta', 'fasta2', 'fasta3']:
            print (os.path.join(outdir_path, chrom, extension))
            os.makedirs(os.path.join(outdir_path, chrom, extension))
        
except Exception as e:
    print (e)         
    


global_coo = pd.read_csv(global_coo_file, sep='\t')[['tr_id', 'strand', 'global_coo_50_and_less']]
global_coo['chr'] = [x.split(':')[0] for x in global_coo['global_coo_50_and_less'].tolist()]
d = global_coo.groupby('chr')[['tr_id', 'global_coo_50_and_less', 'strand']].apply(lambda g: g.values.tolist()).to_dict()


print ('bed files making')

for chrom, v in d.items():
    print ('bed files making', chrom)
    for tr in v:
        tr_id = tr[0]
        coo = tr[1].split('+')
        strand = tr[2]
        chroms = [x.split(':')[0] for x in coo]
        
        # apply -1 to the start 
        starts = [str(int(x.split(':')[1].split('-')[0])-1) for x in coo]
        stops = [x.split(':')[1].split('-')[1] for x in coo]
        
        for i, chr_start_stop in enumerate(zip(chroms, starts, stops)):
            chr_ = chr_start_stop[0]
            start_ = chr_start_stop[1]
            stop_ = chr_start_stop[2]
            idx = tr_id + '_' + str(i)
      
            # /home/DATA2/alla/N_TEMINAL_EXT/paper/data50/chrZ/bed/ENSTxxx.bed
            f = open('%s/%s/bed/%s.bed' % (outdir_path, chrom, idx), 'w')
            f.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (chr_, start_, stop_, idx, '1', strand))
            f.close()        
          
###########################################################################        

def pre_process(chrom):
    print ('start...%s' % chrom)
  
    cmd0 = 'wget https://bds.mpi-cbg.de/hillerlab/120MammalAlignment/Human120way/data/maf/%s.maf.gz -P %s' % (chrom, raw_maf_dir)
    print (chrom, 'wget')
    subprocess.call(cmd0, shell=True)
    
    cmd1 = 'gunzip -c %s/%s.maf.gz > %s/%s.maf' % (raw_maf_dir, chrom, raw_maf_dir, chrom)
    print (chrom, 'gunzip')
    subprocess.call(cmd1, shell=True)
    
    # bb file
    cmd2 = '%s/mafIndex -chromSizes=%s  %s/%s.maf %s/%s.bb' % (mafIndex_dir, chromSizes, raw_maf_dir, chrom, raw_maf_dir, chrom)
    print (chrom, 'bb')
    subprocess.call(cmd2, shell=True)
    
    print (chrom, 'out maf start!')
    for bed_path in glob.glob('%s/%s/bed/*' % (outdir_path, chrom)):
        tmp = pd.read_csv(bed_path, header=None, sep='\t')
        tmp = tmp.sort_values(by=[0, 1, 2])
        tmp.to_csv(bed_path, sep='\t', index=False, header=None)
        cmd3 = '%s/mafExtract %s/%s.bb -outDir %s/%s/maf -regionList=%s' % (mafIndex_dir, raw_maf_dir, 
                                                                            chrom, outdir_path, 
                                                                            chrom, bed_path)
        subprocess.call(cmd3, shell=True)

    
    cmd4 = 'rm -rf %s/%s.maf' % (raw_maf_dir, chrom)
    cmd5 = 'rm -rf %s/%s.maf.gz*' % (raw_maf_dir, chrom)
    cmd6 = 'rm -rf %s/%s.bb*' % (raw_maf_dir, chrom)
    
    print ('rm .maf file')
    subprocess.call(cmd4, shell=True)
    print ('rm .maf.gz file')
    subprocess.call(cmd5, shell=True)
    print ('rm .bb file')
    subprocess.call(cmd6, shell=True)
    

chr_li = ['chrY', 'chrX', 'chr1', 'chr2', 'chr3', 'chr4', 
          'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 
          'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 
          'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']

print ('download maf and extract sub mafs')

if __name__ == '__main__':
    p = Pool(24)
    print(p.map(pre_process, chr_li))
