# Genetic Algorithm for Structural Alignment Cleaning

## Introduction
In the field of homology gene/protein prediction, a crucial task is to identify candidate molecules in target genomes/sequences. This code aims to construct structural alignments from homologous sequences, which can serve as a reference to build predictive models like covariance models. Traditionally, this task involves manual curation, as seen in the Rfam database's alignment of multiple ncRNA families, such as [Let-7](https://rfam.xfam.org/family/RF00027/alignment?acc=RF00027&format=stockholm&download=0).

## Idea
This code utilizes sequences from the miRBase database to assemble multiple structural alignments. It achieves this by selecting sequences based on three parameters:
- Species clade
- Alignment identity
- Sequence quality, determined by the miRBase database

Using these parameters as individuals, a genetic algorithm iterates over the solution space and identifies local maxima based on a fitness score function. After 20 iterations, the output will be the structural alignment with the highest score. The algorithm employs mutation, selection, and crossover functions to generate children, which form the population after each iteration.

## Example: mir-30 Family (MIPF0000005)

The miRBase database initially contains 178 [sequences](https://www.mirbase.org/summary.shtml?fam=MIPF0000005). The objective is to identify sequence groups suitable for constructing structural alignments that optimize the following features:
- Diverse range of species
- Larger alignments
- Lower folding energy
- Correct secondary folding, avoiding patterns like `"\)\.*\("` within the stem region of the alignment.

## Installation and Environment Activation
To run the code, please follow these steps:

1. Clone the repository and navigate to the project directory:
   ```
   git clone git@github.com:cavelandiah/GAlignment.git
   cd GAlignment/
   ```

2. Create the Conda environment using the provided YAML file:
   ```
   conda env create -f geneticalgo.yml
   ```

3. Activate the created environment:
   ```
   conda activate galignment
   ```

## Running the Code
Execute the following commands in the terminal:

1. Create the necessary directories:
   ```
   mkdir Test
   mkdir LOG
   ```

2. Run the code for the desired family (e.g., MIPF0000005):
   ```
   ./evolutionModel_V2.py MIPF0000005 Test LOG
   ```

During the first run, the code may take longer due to the update of the NCBITaxa database. The output will be printed on the screen, showing the results for each iteration. The last line displays the best candidate:

```
#iteration best_score best_individual fitness_detail
1 838.9000000000001 normal,Mammalia,70 [838.9000000000001, 838.9000000000001, 838.9000000000001, 838.9000000000001, 838.9000000000001]
2 838.9000000000001 normal,Mammalia,70 [838.9000000000001, -1000.0, 838.9000000000001, 838.9000000000001, 838.9000000000001]
3 2597.5 high,Mammalia,50 [841.32, 838.9000000000001, 838.9000000000001, 2597.5, 838.9000000000001]
4 2597.5 high,Mammalia,50 [2597.5, 2597.5, 838.9000000000001, 838.9000000000001, 2597.5]
...
20 5954.35 normal,Vertebrata,50 [5954.35, 5954.35, 4249.070000000001, 4244.070000000001, 4244.070000000001]
```

This output indicates that the "GAlignment" algorithm found the best structural alignment by selecting sequences from the mir-30 family, belonging to Vertebrates species, and restricting the alignment to sequences with at least 50% identity and 40% family-sequences. The chosen sequences are from the "normal" category, excluding those from the "high confidence" set.

The resulting alignment can be found in the corresponding file:
```
less -S MIPF0000005_normal_Vertebrata_50.stk
```

The fitness behavior throughout the run can be visualized in the following graph:

![Fitness Graph](https://github.com/cavelandiah/GAlignment/blob/main/MIPF0000005.png)

## Code Length
The code consists of 394 lines, excluding empty lines and comments.
