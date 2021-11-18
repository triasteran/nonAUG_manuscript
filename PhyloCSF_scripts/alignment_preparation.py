import glob
import os
import argparse 
import errno

import subprocess
from subprocess import call
from multiprocessing import Pool

import pandas as pd
import numpy as np
import collections
from collections import defaultdict, OrderedDict

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO,  AlignIO

parser = argparse.ArgumentParser()
parser.add_argument('-a', '--animals_file', required=True) # path to 120way_animals.txt 
parser.add_argument('-g', '--global_coo_file', required=True) # path to global coo file
parser.add_argument('-m', '--main_dir', required=True) # path to main directory 
parser.add_argument('-d', '--dirname', required=True) # name of inner directory
parser.add_argument('-msa', '--msa_view_path', required=True) # path to mas_view tool 
    
args = parser.parse_args()
global_coo_file = args.global_coo_file
animals_file = args.animals_file
main_dir = args.main_dir
dirname = args.dirname 
msa_view_path = args.msa_view_path
outdir_path = os.path.join(main_dir, dirname) 

animals = pd.read_csv(animals_file, sep='\t', header=None)
animals_120 = animals[0].tolist()

animals_120_d = dict(zip(animals[0].tolist(), animals[1].tolist()))

global_coo = pd.read_csv(global_coo_file, sep='\t')[['tr_id', 'strand', 'global_coo_50_and_less']]
global_coo.columns = ['tr_id', 'strand', 'global_coo']
global_coo['chr'] = [x.split(':')[0] for x in global_coo['global_coo'].tolist()]


def merge_subfiles_into_1(path_to_chr_maf, path_out):
    d = defaultdict(list)
    for path in glob.glob(path_to_chr_maf+'*'):
        tr_id = '_'.join(path.split('/')[-1].split('.maf')[0].split('_')[0:-1])   
        if tr_id in d.keys():
            d[tr_id].append(path)
        else:
            d[tr_id] = [path]
        
    for k, v in d.items():
        v1 = sorted(v, key = lambda x: int(x.split('/')[-1].split('.maf')[0].split('_')[-1]))
        new_name = path_out +'/' + '_'.join(v1[0].split('/')[-1].split('.maf')[0].split('_')[0:-1]) + '.maf'
        print (new_name)
        cmd = 'cat %s > %s' % (' '.join(v1), new_name)
        subprocess.call(cmd, shell=True)
        
def stitch_blocks(input_maf_path, output_fasta_path, strand):
    '''
    per 1 transcript read block of alignment
    '''
    blocks =  {}
    all_species_in_maf = set()
    
    for i, multiple_alignment in enumerate(AlignIO.parse(input_maf_path, "maf")):
        blocks[i] = {}
        for seqrec in multiple_alignment:
            blocks[i][seqrec.id] = str(seqrec.seq)
            all_species_in_maf.add(seqrec.id)#
        
    all_species_in_maf = list(all_species_in_maf)
    
    '''
    corrected dict with new blocks 
    '''
    number_of_blocks = len(blocks.keys())
    new_blocks_part = dict.fromkeys(all_species_in_maf , '')
    new_blocks = {i : new_blocks_part for i in list(range(0, number_of_blocks))}
    
    new_blocks1 = { i : {} for i in list(range(0,number_of_blocks))}

    for k, v in new_blocks.items(): #k = 0, 1, 2 - number of block
        len_block = len(blocks[k][list(blocks[k].keys())[0]])
        for k1, v1 in v.items():
            if k1 in blocks[k].keys(): #k1 = rec_id, animal name
                new_blocks1[k][k1] = blocks[k][k1]
            else:
                new_blocks1[k][k1] =  ''.join([str(i) for i in np.repeat('-', len_block)])
                
    #strand + => 0, 1, 2, 3...
    d = {}
    
    if strand == '+':
        for i in range(number_of_blocks): 
            block_ = new_blocks1[i]
            for k1, v1 in block_.items():  
                new_k1 = ''.join(list(filter(lambda x: x.isalpha(), k1.split('.')[0])))
                if new_k1 in d.keys():
                    d[new_k1] += v1 
                else:
                    d[new_k1] = v1
                    
    if strand == '-':
        for i in range(number_of_blocks): 
            block_ = new_blocks1[i]
            for k1, v1 in block_.items():  
                new_k1 = ''.join(list(filter(lambda x: x.isalpha(), k1.split('.')[0])))
                if new_k1 in d.keys():
                    d[new_k1] += v1 
                else:
                    d[new_k1] = v1
        
        # reverse complement seq
        for sp, seq in d.items():
            d[sp] = str(Seq(seq, generic_dna).reverse_complement())
               
    #sort dict: hg38 at the start
    index_map = {v: i for i, v in enumerate(animals_120)}
    
    #return sorted_d
    sorted_d = sorted(d.items(), key=lambda pair: index_map[pair[0]])
    
    #write to file 
    outf = open(output_fasta_path, 'w')
    for el in sorted_d:
        k = el[0]
        v = el[1]
        outf.write('>'+k+'\n')
        outf.write(v+'\n')
    outf.close()
    
def remove_gaps_and_reformat_fasta(input_path, output_path, output_path1):
    cmd = '%s/msa_view  %s --gap-strip ALL > %s' % (msa_view_path, input_path, output_path)
    subprocess.call(cmd, shell=True)
    
    outf = open(output_path1, 'w')
    with open(output_path, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            rec = str(record.id)
            seq = str(record.seq)
            outf.write('>'+animals_120_d[rec]+'\n')
            outf.write(seq+'\n')        
    outf.close()
    
    
print ('merge')
chr_li = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 
          'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 
          'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

for chrom in chr_li:
    print ('merge: ', chrom)
    path_to_chr_maf = os.path.join(outdir_path, chrom, 'maf/')
    path_out = os.path.join(outdir_path, chrom, 'maf_m/')
    print ('path in', path_to_chr_maf)
    print ('path out', path_out)
    merge_subfiles_into_1(path_to_chr_maf=path_to_chr_maf,
                     path_out = path_out)
    
print ('stitch blocks')
for chrom in chr_li:
    print ('stitch: ', chrom)
    names = global_coo[global_coo['chr'] == chrom].tr_id.tolist()
    for name in names:
        strand = global_coo[global_coo['tr_id'] == name].iloc[0].strand
        input_maf_path = os.path.join(outdir_path, chrom, 'maf_m/', '%s.maf' % name)
        output_fasta_path = os.path.join(outdir_path, chrom, 'fasta/', '%s.fasta' % name)
        stitch_blocks(input_maf_path = input_maf_path, 
              output_fasta_path = output_fasta_path, 
              strand = strand)        
        
        
print ('remove gaps')        
for chrom in chr_li:
    print ('remove gaps: ', chrom)
    names = global_coo[global_coo['chr'] == chrom].tr_id.tolist()
    for name in names:
        input_path = os.path.join(outdir_path, chrom, 'fasta/', '%s.fasta' % name)
        output_path = os.path.join(outdir_path, chrom, 'fasta2/', '%s.fasta' % name)
        output_path1 = os.path.join(outdir_path, chrom, 'fasta3/', '%s.fasta' % name)
        remove_gaps_and_reformat_fasta(input_path = input_path, 
                               output_path = output_path, 
                               output_path1 = output_path1)
