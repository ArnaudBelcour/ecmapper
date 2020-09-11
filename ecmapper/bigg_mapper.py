import csv
import re

from Bio import SeqIO

def map_ec_to_bigg(expazy_enzyme_file, bigg_reactions_file, genbank_file):
    # From expazy enzyme, get the name of each EC number.
    # From file: ftp://ftp.expasy.org/databases/enzyme/enzyme.dat
    data_enzymes = {}
    with open(expazy_enzyme_file, 'r') as f:
        enzymes = f.read().split('//')[2:]
        for enzyme_block in enzymes:
            id_enzyme = None
            id_description = None
            for enzyme_line in enzyme_block.split('\n'):
                if 'ID   ' in enzyme_line:
                    id_enzyme = enzyme_line.split('   ')[1]
                if 'DE   ' in enzyme_line:
                    id_description = enzyme_line.split('   ')[1].strip('.')
            if id_enzyme and id_description:
                data_enzymes[id_enzyme] = id_description

    # Map EC name with bigg reactions
    # From file: http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt
    reg = r'(?P<reactants>[a-zA-Z\(\)\d\s\[\]\+\_]*)([\<]*\-[\>]*)(?P<products>[a-zA-Z\(\)\d\s\[\]\+\_]*)'
    reg_expr = re.compile(reg)

    reactions_metabolites = {}
    reactions_ecs = {}
    with open(bigg_reactions_file, 'r') as bigg_reactions_file:
        csvreader = csv.reader(bigg_reactions_file, delimiter='\t')
        next(csvreader)
        for line in csvreader:
            # Extract compounds
            reactants = []
            products = []
            search_result = reg_expr.search(line[2])
            input_reactants = search_result.group('reactants')
            input_products = search_result.group('products')
            for reactant in input_reactants.split(' + '):
                reactant = reactant.strip(' ')
                if len(reactant.split(' ')) >1:
                    reactant = reactant.split(' ')[1]
                reactants.append(reactant)
            for product in input_products.split(' + '):
                product = product.strip(' ')
                if len(product.split(' ')) >1:
                    product = product.split(' ')[1]
                products.append(product)          
            reactions_metabolites[line[0]] = (reactants, products)
            # Extract reactions
            for id_enzyme in data_enzymes:
                if data_enzymes[id_enzyme] in line[1]:
                    if id_enzyme not in reactions_ecs:
                        reactions_ecs[id_enzyme] = [line[0]]
                    else:
                        reactions_ecs[id_enzyme].append(line[0])


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
    reaction_genes = {}
    for gene in metabolic_network:
        for reaction in metabolic_network[gene]:
            if reaction not in reaction_genes:
                reaction_genes[reaction] = [gene]
            else:
                if gene not in reaction_genes[reaction]:
                    reaction_genes[reaction].append(gene)

    return reactions_metabolites, reaction_genes, metabolic_network
