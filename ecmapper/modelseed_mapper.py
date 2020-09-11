import csv
import re

from Bio import SeqIO

def map_ec_to_modelseed(modelseed_reactions_file, genbank_file):
    # Map EC number with modelseed reactions
    # From file: https://github.com/ModelSEED/ModelSEEDDatabase/raw/master/Biochemistry/reactions.tsv
    reg = r'(?P<reactants>[a-z\(\)\d\s\[\]\+]*)([\<]*\=[\>]*)(?P<products>[a-z\(\)\d\s\[\]\+]*)'
    reg_product = r'\(\d\)(?P<compound>\scpd[\d]{5}\[\d\])[\s]*'
    reg_expr = re.compile(reg)
    compound_exp = re.compile(reg_product)
    reactions_ecs = {}
    reaction_compounds = {}
    with open(modelseed_reactions_file, 'r') as bigg_reactions_file:
        csvreader = csv.reader(bigg_reactions_file, delimiter='\t')
        next(csvreader)
        for line in csvreader:
            search_result = reg_expr.search(line[6])
            reactants = []
            for element in search_result.group("reactants").split('+'):
                compound_result = compound_exp.search(element)
                if compound_result is not None:
                    reactants.append(compound_result.group("compound").split('[')[0].strip())
            products = []
            for element in search_result.group("products").split('+'):
                compound_result = compound_exp.search(element)
                if compound_result is not None:
                    products.append(compound_result.group("compound").split('[')[0].strip())
            reaction_compounds[line[0]] = (reactants, products)
            for ec_number in line[13].split('|'):
                if ec_number not in reactions_ecs:
                    reactions_ecs[ec_number] = [line[0]]
                else:
                    reactions_ecs[ec_number].append(line[0])


    # Link genes associated to EC to bigg reactions.
    metabolic_network = {}

    for record in SeqIO.parse(genbank_file, 'genbank'):
        for feature in record.features:
            if feature.type =='CDS':
                if "EC_number" in feature.qualifiers:
                    for ec in feature.qualifiers["EC_number"]:
                        if ec in reactions_ecs:
                            if feature.qualifiers['locus_tag'][0] not in metabolic_network:
                                    metabolic_network[feature.qualifiers['locus_tag'][0]] = reactions_ecs[ec]
                            if feature.qualifiers['locus_tag'][0] not in metabolic_network:
                                metabolic_network[feature.qualifiers['locus_tag'][0]].extend(reactions_ecs[ec])

    # Create a dictionary reaction -> genes.
    reaction_genes={}
    for gene in metabolic_network:
        for reaction in metabolic_network[gene]:
            if reaction not in reaction_genes:
                reaction_genes[reaction] = [gene]
            else:
                if gene not in reaction_genes[reaction]:
                    reaction_genes[reaction].append(gene)

    return reaction_compounds, reaction_genes, metabolic_network
