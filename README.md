# Genetic Algorithm to clean structural alignments
An important task in homology gene/protein prediction is the identification of candidate molecules in target genomes/sequences. In this case I approached the problem in the construction of structural alignments from homologous sequences. This task is performed usually by a hand curation, see for example the Rfam database [Let-7](https://rfam.xfam.org/family/RF00027/alignment?acc=RF00027&format=stockholm&download=0). 

This code uses the sequences from the miRBase database, and assembles multiple structural alignments by playing with the following variables:
Target clade, Identity of the alignment and quality of the sequences. The final output will be the structural alignment that scores a maximum after 20 iterations. 

## Example: Family mir-30 (MIPF0000005)

By default the miRBase database has 178 [sequences](https://www.mirbase.org/summary.shtml?fam=MIPF0000005). The idea is to identify the group of sequences that are suitable to build the structural alignment.

Please run as:

    conda install -f geneticalgo.yml
    git clone git@github.com:cavelandiah/GAlignment.git
    cd GAlignment/
    mkdir Test
    mkdir LOG
    ./evolutionModel_V2.py MIPF0000005 Test LOG
    
The output will printed in the screen, indicating at the last line the best structural alignment.
