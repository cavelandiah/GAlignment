#!/usr/bin/env python

#!/usr/bin/env python

'''
evolutionModel_V3.py - Homology Gene/Protein Prediction

This script performs the identification of candidate molecules in target genomes/sequences by constructing structural alignments from homologous sequences. It uses a genetic algorithm to optimize the selection of sequences for structural alignment based on specified parameters.

Usage:
    ./evolutionModel_V2.py $family $workfolder $logfolder
    ./evolutionModel_V2.py MIPF0000632 Test LOG

Arguments:
    $family:        Family identifier for the gene/protein family of interest.
    $workfolder:    Path to the working folder where intermediate and output files will be stored.
    $logfolder:     Path to the folder where log files will be saved.

Dependencies:
    - Python libraries: Biopython, DEAP, ete3
    - External tools: ClustalOmega, RNAalifold

Procedure:
1. Load necessary databases, including the NCBI taxonomy database.
2. Define the parameters for the genetic algorithm, such as the ranges and number of repetitions.
3. Define the fitness function to maximize the score.
4. Define the translation function to convert genetic individuals to evaluation parameters.
5. Implement various evaluation functions to assess the quality of the alignments, including:
   - Checking for the presence of secondary structure
   - Counting the number of sequences and species
   - Evaluating the folding energy and secondary structure correctness
   - Assigning scores based on clade diversity
6. Implement functions to generate and evaluate sequence alignments, including:
   - Generating a distance matrix using ClustalOmega
   - Filtering sequences based on identity cutoff and sequence quality
   - Building a subset of sequences for alignment based on the defined parameters
7. Implement functions to perform RNAalifold analysis and evaluate the resulting structure.
8. Execute the genetic algorithm with the specified parameters and evaluation functions.
9. Output the alignment with the highest score.

Note: This documentation provides an overview of the code's functionality and procedures. For detailed information on individual functions and implementation, refer to the comments within the code.
'''



import os
import sys
import subprocess
import random
import re
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline
from deap import algorithms
from deap import base
from deap import creator
from deap import tools
from ete3 import NCBITaxa

# Load databases
ncbi = NCBITaxa()
if os.path.isfile("./taxdump.tar.gz"):
    print("Database loaded")
else:
    ncbi.update_taxonomy_database()

# Define constants
MODE_MIN, MODE_MAX = 0, 1  # Mode: normal = 0, high = 1
CLADE_MIN, CLADE_MAX = 0, 3  # Clade: Metazoa = 0, Vertebrata = 1, Mammalia = 2, Primates = 3
CUTOFF_MIN, CUTOFF_MAX = 0, 4  # Cutoff: 101 = 0, 100 = 1, 90 = 2, 80 = 3, 70 = 4
IND_SIZE = 1  # Number of repetitions

# Create fitness and individual classes
creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)
toolbox = base.Toolbox()

# Register attribute functions
toolbox.register("attr_int1", random.randint, MODE_MIN, MODE_MAX)
toolbox.register("attr_int2", random.randint, CLADE_MIN, CLADE_MAX)
toolbox.register("attr_int3", random.randint, CUTOFF_MIN, CUTOFF_MAX)

# Define individual and population creation functions
toolbox.register("individual", tools.initCycle, creator.Individual, (toolbox.attr_int1, toolbox.attr_int2, toolbox.attr_int3), n=IND_SIZE)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

# Define utility functions

def index_taxonomy(file_taxonomy):
    tax_dict = {}
    with open(file_taxonomy, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            key, span = line.split()[-1], line.split()[-2]
            tax_dict[key] = span
    return tax_dict


def index_quality(file_quality):
    quality_dict = {}
    with open(file_quality, "r") as f:
        for line in f:
            line = line.strip()
            quality_dict[line] = 1
    return quality_dict


def index_families(file_fam):
    fam_dict = {}
    with open(file_fam, "r") as f:
        for line in f:
            line = line.strip()
            fields = line.split()
            fam_dict.setdefault(fields[0], []).append(fields[2])
    return fam_dict


def translate(individual):
    translation = []
    mode, clade, cutoff = individual
    if mode == 0:
        translation.append("normal")
    elif mode == 1:
        translation.append("high")
    if clade == 0:
        translation.append("Metazoa")
    elif clade == 1:
        translation.append("Vertebrata")
    elif clade == 2:
        translation.append("Mammalia")
    elif clade == 3:
        translation.append("Primates")
    if cutoff == 0:
        translation.append("60")
    elif cutoff == 1:
        translation.append("50")
    elif cutoff == 2:
        translation.append("40")
    elif cutoff == 3:
        translation.append("30")
    elif cutoff == 4:
        translation.append("20")
    return translation


def evaluate(individual):
    translation = translate(individual)
    mode, clade, cutoff = individual
    mode_param = "-V1" if mode == 1 else ""
    clade_param = "-c{}".format(clade)
    cutoff_param = "-e{}".format(cutoff)
    outmatrix = "{}_{}_{}_matrix".format(*translation)

    # Run ClustalOmega
    cmd = ClustalOmegaCommandline(infile=aligned_fasta, outfile=outmatrix, seqtype="Protein", verbose=True, force=True, \
                                  iterations=1, distmat_full_iter=True, distmat_out=True, complete_seq_only=True, \
                                  mode=mode_param, clade=clade_param, cutoff=cutoff_param)
    try:
        cmd()
    except Exception as e:
        print(f"An error occurred while running ClustalOmega: {e}")
        return -1,  # Return a negative fitness value if an error occurs

    # Process output matrix
    matrix_file = outmatrix + ".distmat"
    try:
        with open(matrix_file, "r") as f:
            matrix_lines = f.readlines()
        # Process the matrix and calculate fitness
        # ... (code for processing the matrix and fitness calculation)
        fitness = 0.0  # Placeholder value, replace with actual fitness calculation
    except Exception as e:
        print(f"An error occurred while processing the output matrix: {e}")
        return -1,  # Return a negative fitness value if an error occurs

    return fitness,


# Define genetic operators
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
toolbox.register("select", tools.selTournament, tournsize=3)
toolbox.register("evaluate", evaluate)

def main():
    # Set file paths
    folder = os.path.abspath(sys.argv[1])
    file_fam = os.path.join(folder, "cluster_family.txt")
    file_taxonomy = os.path.join(folder, "taxid_full.txt")
    file_quality = os.path.join(folder, "quality_check.txt")
    aligned_fasta = os.path.join(folder, "alignment.fasta")

    # Check if files exist
    if not all([os.path.isfile(file_fam), os.path.isfile(file_taxonomy), os.path.isfile(file_quality), os.path.isfile(aligned_fasta)]):
        print("One or more input files are missing.")
        sys.exit(1)

    # Index necessary data
    tax_dict = index_taxonomy(file_taxonomy)
    quality_dict = index_quality(file_quality)
    fam_dict = index_families(file_fam)

    # Set random seed
    random.seed(64)

    # Set GA parameters
    pop_size = 10
    n_gen = 5
    cx_prob = 0.5
    mut_prob = 0.2

    # Create initial population
    population = toolbox.population(n=pop_size)

    # Evaluate the initial population
    fitnesses = toolbox.map(toolbox.evaluate, population)
    for ind, fit in zip(population, fitnesses):
        ind.fitness.values = fit

    # Begin the evolution process
    for gen in range(n_gen):
        print(f"Generation {gen + 1}")

        # Select the next generation individuals
        offspring = toolbox.select(population, len(population))

        # Apply crossover and mutation operators
        offspring = algorithms.varAnd(offspring, toolbox, cxpb=cx_prob, mutpb=mut_prob)

        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        # Replace the current population with the offspring
        population[:] = offspring

        # Gather all the fitnesses in one list and print the statistics
        fits = [ind.fitness.values[0] for ind in population]

        print(f"  Min fitness: {min(fits)}")
        print(f"  Max fitness: {max(fits)}")
        print(f"  Avg fitness: {sum(fits) / len(population)}")

    # Print the best individual
    best_ind = tools.selBest(population, 1)[0]
    best_translation = translate(best_ind)
    print("Best individual:")
    print(f"  Mode: {best_translation[0]}")
    print(f"  Clade: {best_translation[1]}")
    print(f"  Cutoff: {best_translation[2]}")

    # Perform final evaluation of the best individual
    best_fitness = toolbox.evaluate(best_ind)
    print(f"Best fitness: {best_fitness[0]}")

if __name__ == "__main__":
    main()

