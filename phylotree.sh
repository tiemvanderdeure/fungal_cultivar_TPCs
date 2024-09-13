mafft --ep 0 --genafpair --maxiterate 1000 input sequences.fasta > alignment.fasta
iqtree -s alignment.fasta -m TEST -bb 1000 -alrt 1000
