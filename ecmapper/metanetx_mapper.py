import csv
import re

from Bio import SeqIO


def extract_formula(reaction_formula):
    reg = r'(?P<reactants>[A-Z\-\+a-z\(\)\d\s\[\]\+\@,]*)(\=)(?P<products>[A-Z\-\+a-z\(\)\d\s\[\]\,@]*)'
    reg_product = r'\d\s(?P<compound>[A-Z\d@]*)'
    reg_expr = re.compile(reg)
    compound_exp = re.compile(reg_product)

    search_result = reg_expr.search(reaction_formula)
    reactants = []
    for element in search_result.group("reactants").split(' + '):
        compound_result = compound_exp.search(element)
        reactants.append(compound_result.group("compound"))
    products = []
    for element in search_result.group("products").split(' + '):
        compound_result = compound_exp.search(element)
        if compound_result is not None:
            products.append(compound_result.group("compound"))
    return reactants, products


def map_ec_to_metanetx(metanetx_reactions_file, genbank_file):
    # Extract reactions and metabolites from metanetx.
    reactions_metabolites = {}
    reactions_ecs = {}
    ex_reaction_ecs = {}
    with open(metanetx_reactions_file, 'r') as bigg_reactions_file:
        csvreader = csv.reader(bigg_reactions_file, delimiter='\t')
        next(csvreader)
        for line in csvreader:
            if not line[0].startswith('EMPTY') and not line[0].startswith('#'):
                if line[3] != '':
                    for ec in line[3].split(','):
                        if ec not in reactions_ecs:
                            reactions_ecs[ec] = [line[0]]
                        else:
                            reactions_ecs[ec].append(line[0])
                        if ec not in ex_reaction_ecs:
                            ex_reaction_ecs[ec] = [line[2]]
                        else:
                            ex_reaction_ecs[ec].append(line[2])
                reactions_metabolites[line[0]] = extract_formula(line[1])


    # Extract genes from genbank file.
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

    return reactions_metabolites, reaction_genes, metabolic_network
