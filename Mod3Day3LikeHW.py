# GOAL: from a population of solutions for a given set of loci, score the genes on the loci 
# using the method in Tasan et al. Once scored, visualize the gene scores and subnetwork 
# and induce on the set of loci in Figure 1 of the Tasan et al. paper. 

import pandas as pd
import numpy as np
import random

string = pd.read_csv(file_name, delimiter='\t', names=['Gene1', 'Gene2', 'Weight'])

# define a function to create dictionary with locus as keys and genes as values
def fa_genes(file_name):
    FA_genes = {}
    for line in open(file_name, 'r'):
        locus = line.split("\t")[0] 
        genes = line.split("\t")[2:]
        FA_genes[locus] = genes 
    return FA_genes

FA_genes = fa_genes(file_name)

# define a function to create a network of FA genes
def FA_network(FA_genes, string, Gene1, Gene2):
    FA_genes = set(gene for genes in FA_genes.values() for gene in genes)
    mask = string[Gene1].isin(FA_genes) & string[Gene2].isin(FA_genes)
    network = string.loc[mask, [Gene1, Gene2]].values.tolist()
    return network

FA_network = FA_network(FA_genes, string, 'Gene1', 'Gene2')

# convert FA_network to a dataframe
FA_df = pd.DataFrame(FA_network, columns=['Gene1', 'Gene2'])

# random FA subnetwork with 1 gene per locus
def random_subnetworks(FA_genes):
    random_nets = {locus: random.choice(FA_genes[locus]) for locus in FA_genes.keys()}
    return [[random.choice(FA_genes[locus]) for locus in FA_genes.keys()] for _ in range(5000)]

prix_fixe = random_subnetworks(FA_genes)

def gene_scoring(prix_fixe, FA_genes, num_networks):
    def compute_density(subnetwork, full_network):
        subnet_df = full_network.loc[(full_network['Gene1'].isin(subnetwork)) & (full_network['Gene2'].isin(subnetwork))]
        return len(subnet_df)

    gene_scores = {}
    loci = list(FA_genes.keys())
    full_densities = [compute_density(network, FA_df) for network in prix_fixe]
    prix_fixe = [set(network) for network in prix_fixe]  # Convert to sets for faster membership tests

    for i in range(num_networks):
        if len(prix_fixe[i]) < len(loci):
            print(f"Network {i} has fewer genes than expected.")
            continue
        for j, locus in enumerate(loci):
            g_star = random.choice(list(prix_fixe[i]))  # Select a random gene from the subnetwork
            new_network = prix_fixe[i].copy() # make a new network with g_star removed
            new_network.remove(g_star) # remove g_star from the new network
            empty_density = compute_density(new_network, FA_df) # calculate the density of the new network
            density_gi = full_densities[i] - empty_density  # calculate contribution of g* to the connectivity of the loci
            gene_scores[g_star] = density_gi / len(loci) # record contribution of g* to the connectivity of the loci
    return gene_scores

gene_scores = gene_scoring(prix_fixe, FA_genes, 5000)

# convert gene_scores to a dataframe
gene_scores_df = pd.DataFrame(gene_scores.items(), columns=['Gene', 'Score'])

# Create a reverse lookup dictionary
reverse_dict = {gene: key for key, value in FA_genes.items() for gene in value}

# Use the reverse lookup dictionary to append the locus
gene_scores_df['Locus'] = gene_scores_df['Gene'].map(reverse_dict)

# export FA_df to a txt file
FA_df.to_csv(file_name, sep='\t', index=False)

#export gene_scores_df to a txt file
gene_scores_df.to_csv(file_name, sep='\t', index=False)

