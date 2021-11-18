function parallel_call {
    dir1=data # path to directory with data 
    PH=PhyloCSF # path to PhyloCSF software 
    echo $1;
    cd  $dir1/$1/fasta3;
    parallel --jobs 10 --keep-order $PH/PhyloCSF 120mammals --removeRefGaps \
     $dir1/$1/fasta3/{} > \
     $dir1/$1/PhyloCSF_out.txt ::: *.fasta
}

export -f parallel_call

chrom_set='chr20 chr21 chr22 chrX chrY chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10'
parallel --verbose -j 3 parallel_call ::: ${chrom_set}
