# Genetic Algorithm to clean structural alignments
## Context
An important task in homology gene/protein prediction is the identification of candidate molecules in target genomes/sequences. I approached the construction of structural alignments from homologous sequences, which would serve as reference to build predictive models as covariance models. 
This task is performed usually by a hand curation, as seen in the Rfam database [Let-7](https://rfam.xfam.org/family/RF00027/alignment?acc=RF00027&format=stockholm&download=0) for multiple ncRNA families.  

## Idea
This code uses the sequences from the miRBase database, and assembles multiple structural alignments by selecting sequences to build structural alignments using three parameters:
- Species clade
- Identity of the alignment 
- Quality of the sequences, based on the miRBase database

Based on the codification of those parameters as individuals, a genetic algorithm offered the way to iterate over the iteration space and locate local maxima, based on a score function (fitness score). The final output will be the structural alignment that scores a maximum after 20 iterations. It uses mutation, selection and cross functions to generate the children that would be the population after each iteration.  

## Example: Family mir-30 (MIPF0000005)

By default the miRBase database has 178 [sequences](https://www.mirbase.org/summary.shtml?fam=MIPF0000005). The idea is to identify groups of sequences that are suitable to build the structural alignments that optimize the following features:

- Diverse number of species
- Larger alignments
- Lower folding energy
- Correct secondary folding, avoiding patters as `"\)\.*\("` inside the stem region of the alignment.

## Instalation and environment activation
To do so, please run as:

    git clone git@github.com:cavelandiah/GAlignment.git
    cd GAlignment/
    conda install -f geneticalgo.yml
    conda activate galignment
### Running 
    mkdir Test
    mkdir LOG
    ./evolutionModel_V2.py MIPF0000005 Test LOG
    
When ran for the first time it takes longer due update of the NCBITaxa database. The output will printed in the screen, indicating the results for each iteration. At the last line the best candidate:

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

<img src="https://github.com/cavelandiah/GAlignment/blob/main/MIPF0000005.png" alt="drawing" width="350"/>

### Code length
    
    cat evolutionModel_V2.py | grep -v "^$" | grep -vc "^\#"
    394
    
