# Genetic Algorithm to clean structural alignments
An important task in homology gene/protein prediction is the identification of candidate molecules in target genomes/sequences. In this case I approached the problem in the construction of structural alignments from homologous sequences. This task is performed usually by a hand curation, see for example the Rfam database [Let-7](https://rfam.xfam.org/family/RF00027/alignment?acc=RF00027&format=stockholm&download=0). 

This code uses the sequences from the miRBase database, and assembles multiple structural alignments by playing with the following variables:
Target clade, Identity of the alignment and quality of the sequences. The final output will be the structural alignment that scores a maximum after 20 iterations. 

## Example: Family mir-30 (MIPF0000005)

By default the miRBase database has 178 [sequences](https://www.mirbase.org/summary.shtml?fam=MIPF0000005). The idea is to identify the group of sequences that are suitable to build the structural alignment.

Please run as:

    git clone git@github.com:cavelandiah/GAlignment.git
    cd GAlignment/
    conda install -f geneticalgo.yml
    conda activate galignment
    mkdir Test
    mkdir LOG
    ./evolutionModel_V2.py MIPF0000005 Test LOG
    
The output will printed in the screen, indicating the results for each iteration. At the last line the best candidate:

    #iteration best_score best_individual fitness_detail
    1 838.9000000000001 normal,Mammalia,70 [838.9000000000001, 838.9000000000001, 838.9000000000001, 838.9000000000001, 838.9000000000001]
    2 838.9000000000001 normal,Mammalia,70 [838.9000000000001, -1000.0, 838.9000000000001, 838.9000000000001, 838.9000000000001]
    3 2597.5 high,Mammalia,50 [841.32, 838.9000000000001, 838.9000000000001, 2597.5, 838.9000000000001]
    4 2597.5 high,Mammalia,50 [2597.5, 2597.5, 838.9000000000001, 838.9000000000001, 2597.5]
    ...
    20 5954.35 normal,Vertebrata,50 [5954.35, 5954.35, 4249.070000000001, 4244.070000000001, 4244.070000000001]

It means that ``GAlignment`` found that taking the sequences from the family mir-30 which belongs from Vertebrates species, restricting the alignment for sequences showing 50% identity with at least 40% family-sequences and selecting those sequences from the _normal_ and not from the _high confidence_ set.

The result is in the correspondent file:
    
    less -S MIPF0000005_normal_Vertebrata_50.stk

In terms of fitness, the behaviour throught this run can be visualized as follows:

![MIPF0000005](https://github.com/cavelandiah/GAlignment/blob/main/MIPF0000005.png)
