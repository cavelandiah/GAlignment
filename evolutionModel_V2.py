#!/usr/bin/env python

'''
./evolutionModel_V2.py $family $workfolder $logfolder
'''
import os
import os.path
import sys
import subprocess
import random
import re
#import RNA

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

# In this case I'd like to maximise the score (fitness function)
creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)

toolbox = base.Toolbox()

# Variable codification. Meaning and limits:  <10-11-21, cavelandiah> #
MODE_MIN, MODE_MAX = 0, 1  # Mode: normal = 0, high = 1
CLADE_MIN, CLADE_MAX = 0, 3  # Clade = Metazoa =0, Vertebrata=1, Mammalia=2, Primates=3
CUTOFF_MIN, CUTOFF_MAX = 1, 4  # 101=0,100=1,90=2,80=3,70=4
IND_SIZE = 1  # Number of repetitions

# Register variables = alias, function_to_alias
toolbox.register("attr_int1", random.randint, MODE_MIN, MODE_MAX)
toolbox.register("attr_int2", random.randint, CLADE_MIN, CLADE_MAX)
toolbox.register("attr_int3", random.randint, CUTOFF_MIN, CUTOFF_MAX)
# Build the individual as = [ int, int, int]
toolbox.register("individual", tools.initCycle, creator.Individual, (toolbox.attr_int1, toolbox.attr_int2, toolbox.attr_int3), n=IND_SIZE)
# Population
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

## Def
def index_taxonomy(file_taxonomy):
    tax_dict = {}
    for line in open(file_taxonomy, "r"):
        line = line.strip()
        if line.startswith("#"):
            continue
        key = line.split()[-1]  # 6085
        span = line.split()[-2]  # Metazoa;Bilateria;Deuterostoma;Hemichordata;
        tax_dict[key] = span
    return tax_dict


def index_quality(file_quality):
    quality_dict = {}
    for line in open(file_quality, "r"):
        line = line.strip()
        quality_dict[line] = 1
    return quality_dict


def index_families(file_fam):
    fam_dict = dict()
    for line in open(file_fam, "r"):
        line = line.strip()
        fields = line.split()
        fam_dict.setdefault(fields[0], []).append(fields[2])
    return fam_dict


def translate(individual):
    # Translate the meaning of the individual into parameters of evaluation program
    translation = []
    for i in range(len(individual)):
        if i == 0:  # mode
            if individual[i] == 0:
                translation.append("normal")
            elif individual[i] == 1:
                translation.append("high")
        elif i == 1:  #clade
            if individual[i] == 0:
                translation.append("Metazoa")
            elif individual[i] == 1:
                translation.append("Vertebrata")
            elif individual[i] == 2:
                translation.append("Mammalia")
            elif individual[i] == 3:
                translation.append("Primates")
        elif i == 2:  #cutoff
            if individual[i] == 0:
                translation.append("101")
            elif individual[i] == 1:
                translation.append("100")
            elif individual[i] == 2:
                translation.append("90")
            elif individual[i] == 3:
                translation.append("80")
            elif individual[i] == 4:
                translation.append("70")
    return translation


def evaluate_if_no_structure(str):
    # Evaluate if SS_consensus has structure
    str_len = len(str)
    count_gaps = 0
    for i in str:
        match_gap = re.match(r'\.', i)
        if match_gap:
            count_gaps = count_gaps + 1
    if str_len == count_gaps:
        return 1
    return 0


def file_len(fname):
    # Account number of species + number of sequences
    species = dict()
    lines_number = 0
    with open(fname, 'r') as infile:
        for line in infile:
            line = line.strip()
            if re.search(r'^#', line):
                continue
            else:
                seq_header = line.split()[0]
                spe = seq_header.split('-')[0]
                species[spe] = 1
                lines_number = lines_number + 1
    number_species = len(species.keys())
    # score = (1 / int(lines_number) / int(number_species))
    score = number_species + lines_number
    return score


def evaluate_energy(stofile):
    # Evaluate folding energy of the sto file
    first = ["RNAalifold", "--noPS", "-q", "-r", stofile]  # -r = RIBOSUM scoring
    f = subprocess.check_output(first)
    if not isinstance(f, str):
        f = f.decode("utf-8")
    energy = f.split(' (')[-1]
    energy = energy.split('=')[0]
    energy = float(energy)*-1  # Positive to valid energies
    return energy


def evaluate_structure(structure):
    '''
    Evaluate the presence of ...((((...))))..
    penalyze the presence of ..(((())(...))((.))).
    '''
    pattern = "\)\.*\("
    parts_str = len(re.findall(pattern, structure))
    score = int(parts_str) * -10
    return score


def evaluate_clade(clade):
    # Reward some clades in order to have more diversity:  <10-11-21, cavelandiah> #
    score = None
    if clade == "Metazoa":
        score = 10
    elif clade == "Vertebrata":
        score = 5
    elif clade == "Mammalia":
        score = 2.5
    elif clade == "Primates":
        score = 1.25
    else:
        score = 0
    return score


def stockholm_evaluation(stoFile, clade):
    # Evaluate STO alignment:  <10-11-21, cavelandiah> #
    if not os.path.isfile(stoFile):
        return -100
    number_of_align_sequences = file_len(stoFile)
    if number_of_align_sequences == 0:
        return -100
    align = AlignIO.read(stoFile, "stockholm")
    # Obtain information from secondary structure
    structure = align.column_annotations['secondary_structure']  # Obtain SStr
    eval_empty = evaluate_if_no_structure(structure)
    if eval_empty == 1:
        # Here the str is not valid
        return -100
    else:
        # Look number sequences AND energy
        clade_contr = evaluate_clade(clade) # weight based on clade
        folding_energy = evaluate_energy(stoFile)
        mirna_folding = evaluate_structure(structure) # Negative value if found additional stems in str
        # old
        # score = float(number_of_align_sequences) + float(folding_energy) + float(mirna_folding)
        # Here is the fitness score ## EMPIRICAL# See formula:  <10-11-21, cavelandiah> #
        score = float(clade_contr) + (float(number_of_align_sequences)*float(folding_energy)) + float(mirna_folding)
        return score


def build_fasta_file(variables, family, sequences, taxonomy, quality, fam_relation):
    mode = variables[0]
    clade = variables[1]
    cutoff = variables[2]
    seq_dict={}
    fasta_seq = SeqIO.parse(open(sequences), 'fasta')
    ids_family = fam_relation[family]
    #tax_id_selection = ncbi.get_taxid_translator([clade])
    for fasta in fasta_seq:
        vector_selection = [0, 0, 0]
        name, seq = fasta.description, fasta.seq
        header = name.split()
        infer_specie = header[2]+" "+header[3]
        id_mirbase = header[1]
        # Check seq in family
        if id_mirbase in ids_family:
            tax_id = ncbi.get_name_translator([infer_specie])
            tax_id = tax_id.get(infer_specie)[0]
            #lineage = ncbi.get_lineage(tax_id[0])

            # Lineage selection
            if str(tax_id) in taxonomy.keys():
                lineage = taxonomy.get(str(tax_id))
                if clade.lower() in lineage.lower():
                    vector_selection[0] = 1
            # High quality
            if id_mirbase in quality:
                vector_selection[1] = 1
            vector_string = ''.join(str(e) for e in vector_selection)
            seq_dict.setdefault(vector_string, []).append(id_mirbase)

    # Print selected [1, 1/0, 1 ]
    mode_numb = 0
    if mode == "high":
        mode_numb = 1
    # Temporal
    subset = [1, mode_numb, 0]
    string_subset = "".join(str(f) for f in subset)
    selected = seq_dict[string_subset]
    out_name = "selected.fa"
    out_selected = open(out_name,'w')
    fasta_seq2 = SeqIO.parse(open(sequences), 'fasta')
    for fasta in fasta_seq2:
        name, seq = fasta.description, fasta.seq
        header = name.split()
        id_mirbase = header[1]
        if id_mirbase in selected:
            out_selected.write(">"+str(name)+"\n"+str(seq)+"\n")
    out_selected.close()
    return out_name


def evaluate(family, output_folder, output_folder_complete, logfolder, taxonomy, quality, mapping_file, individual):
    # Evaluate the fitness:  <10-11-21, cavelandiah> #
    # ./validate_miRNA_candidates_V2.pl $family $cutoff $clade miRBaseModels_Metazoa_tunicata $mode
    result = None
    values = translate(individual)
    # Individual meaning:  <10-11-21, cavelandiah> #
    # [mode, clade, cutoff]
    errortemp = open(logfolder+"/outerrorNOVALID_"+family+"_"+values[1]+"_"+values[2]+"_"+values[0]+".txt", 'w')
    # outtemp = "./LogsmiRBaseModelsMetazoaJoined/outNOVALID_"+family+"_"+values[1]+"_"+values[2]+"_"+values[0]+".txt"
    #first = ["./validate_miRNA_candidates_V2.pl", family, values[2], values[1], output_folder, values[0], output_folder_complete+"/.."]
    out_file = str(family)+".sto"
    fasta_file_subset = build_fasta_file(values, family, hairpin_fasta, taxonomy, quality, mapping_file)
    # Clustal Omega
    clustalomega_cline = ClustalOmegaCommandline(infile=fasta_file_subset, outfile=out_file, outfmt='st', verbose=False, auto=False)
    clustalomega_cline()
    #sto_file_to_create = output_folder_complete+"/"+family+"_"+values[1]+"_"+values[2]+"_"+values[0]+"/Output/"+family+".out/"+family+".stk"
    if os.path.isfile(out_file):
        result = out_file
    else:
        #folder = output_folder_complete+"/"+family+"_"+values[1]+"_"+values[2]+"_"+values[0]+"/Output/"+family+".out/"
        #if os.path.isdir(folder) and len(os.listdir(folder)) == 0:
        result = 'NA'
        #else:
        #    result = subprocess.check_output(first, stderr=errortemp)
    #score = None
    #if not isinstance(result, str):
    #    result = result.decode("utf-8")
    if re.search(r'NA', result):
        # The evaluation failed, penalize it!
        score = -1000
    else:
        # Get evaluation over sto file in 'result' variable
        score = stockholm_evaluation(result, values[1])
    errortemp.close()
    return score,


def mutVector(individual, indp, MODE_MIN, MODE_MAX, CLADE_MIN, CLADE_MAX, CUTOFF_MIN, CUTOFF_MAX):
    # [a, b, c]
    # a = [0,1]
    # b = [0,1,2,3]
    # c = [0,1,2,3,4]
    for i in range(len(individual)):
        # The mutation threshold is the same as defined
        if random.random() < indp:
            if i == 0:  # Mode
                individual[i] = random.randint(MODE_MIN, MODE_MAX)
            elif i == 1:  # Clade
                individual[i] = random.randint(CLADE_MIN, CLADE_MAX)
            elif i == 2:  # Cutoff
                individual[i] = random.randint(CUTOFF_MIN, CUTOFF_MAX)
    return individual,


def crossVector(child1, child2):
    new_child1 = [child1[0], child2[1], child1[2]]
    new_child2 = [child2[0], child1[1], child2[2]]
    return new_child1, new_child2


def analyse_winners(vector):
    limit = 5
    if len(vector) < limit:
        return 0
    elif len(vector) == limit:
        lista = dict()
        for i in range(len(vector)):
            concatenate = ",".join(vector[i])
            lista[concatenate] = 1
        uniq_set = len(lista.keys())
        if (uniq_set > 1):
            return 0
        elif (uniq_set == 1):
            return 2
        elif (uniq_set < 1):
            return 0
    elif len(vector) > limit:
        vector = vector[-limit:] # restrict to last 5 elements
        lista = dict()
        for i in range(len(vector)):
            concatenate = ",".join(vector[i])
            lista[concatenate] = 1
        uniq_set = len(lista.keys())
        if (uniq_set > 1):
            return 0
        elif (uniq_set == 1):
            return 2
        elif (uniq_set < 1):
            return 0

# Parameters
# toolbox.register("mate", tools.cxUniform)

## Init
family = sys.argv[1]
output_folder_complete = sys.argv[2]
logfolder = sys.argv[3]
output_folder = output_folder_complete.split('/')[-1] #Last part of output folder
taxonomy = index_taxonomy("organisms.txt")
quality = index_quality("miRNA_high_v22.list")
mapping_file = index_families("mapping_file_original_cleaned_mirbase_v22_families.txt")
hairpin_fasta = "hairpin-metazoa.fa"

toolbox.register("mate", crossVector)
toolbox.register("mutate", mutVector)  # , (toolbox.attr_int1, toolbox.attr_int2, toolbox.attr_int3))
toolbox.register("select", tools.selTournament, tournsize=39)
# toolbox.register("select", tools.selRoulette)
toolbox.register("evaluate", evaluate, family, output_folder, output_folder_complete, logfolder, taxonomy, quality, mapping_file)

# Population and start fitness
pop = toolbox.population(n=40)
fitnesses = list(map(toolbox.evaluate, pop))
for ind, fit in zip(pop, fitnesses):
    ind.fitness.values = fit

# Step switches
g=0
switch=0
# Cross child:  <10-11-21, cavelandiah> #
CXBP=0.7
# Mutation in individual:  <10-11-21, cavelandiah> #
MUTPB=0.1

# All fitness for the population:  <10-11-21, cavelandiah> #
fits = [ind.fitness.values[0] for ind in pop]

logs = open(logfolder+"/"+family+"_log.txt", 'a')
selected_winner = []
#while g < 100: #max(fits) < 8:
# Iteration steps, until the same winning is conserved over 20 rounds:  <10-11-21, cavelandiah> #
while switch < 1 and g < 20:
    g = g + 1
    # print("-- Generation %i --" % g)
    # Select the next generation individuals
    offspring = toolbox.select(pop, len(pop))
    # Clone the selected individuals
    offspring = list(map(toolbox.clone, offspring))

#    #Apply crossover on the offspring
    for child1, child2 in zip(offspring[::2], offspring[1::2]):
        if random.random() < CXBP:
            # P(cross both childs < 0.7)
            toolbox.mate(child1, child2)
            del child1.fitness.values
            del child2.fitness.values

    # Apply mutation on the offsprint
    for mutant in offspring:
        if random.random() < MUTPB:
            # The P(mutate_individual < 0.1) and each site 0.5
            toolbox.mutate(mutant, 0.5, MODE_MIN, MODE_MAX, CLADE_MIN, CLADE_MAX, CUTOFF_MIN, CUTOFF_MAX)
            del mutant.fitness.values

    # Evaluate the individuals with an invalid fitness
    invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
    fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
    for ind, fit in zip(invalid_ind, fitnesses):
        ind.fitness.values = fit

    # The population is entirely replaced by the offspring
    pop[:] = offspring

    # Gather all the fitnesses in one list and print the stats
    fits = [ind.fitness.values[0] for ind in pop]
    # STATS of generations:  <10-11-21, cavelandiah> #
    length = len(pop)
    mean = sum(fits) / length
    sum2 = sum(x*x for x in fits)
    std = abs(sum2 / length - mean**2)**0.5
    maximum = max(fits)
    winner = fits.index(maximum)
    winnerTr = translate(pop[winner])
    #print("  Min %s" % min(fits))
    #print("  Max %s" % max(fits))
    #print("  Avg %s" % mean)
    #print("  Std %s" % std)
    #print(fits)
    newline = ",".join(winnerTr)
    logs.write(str(g)+" "+str(newline)+" "+str(maximum)+"\n")
    selected_winner.append(winnerTr)
    switch = analyse_winners(selected_winner)
    if maximum < 0:
        logs.write("No viable alignment")
        sys.exit()
