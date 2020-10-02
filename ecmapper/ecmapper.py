import csv
import os

from ecmapper.bigg_mapper import map_ec_to_bigg
from ecmapper.modelseed_mapper import map_ec_to_modelseed
from ecmapper.metanetx_mapper import map_ec_to_metanetx
from ecmapper.sbml import create_sbml


def workflow(genbank_file, database_folder, output_folder):
    enzyme_file = database_folder + "/enzyme.dat"
    bigg_file = database_folder + "/bigg_models_reactions.txt"
    modelseed_file = database_folder + "/reactions.tsv"
    metanetx_file = database_folder + "/reac_prop.tsv"

    compounds_bigg, reactions_bigg, genes_bigg = map_ec_to_bigg(enzyme_file, bigg_file, genbank_file)

    compounds_modelseed, reactions_modelseed, genes_modelseed = map_ec_to_modelseed(modelseed_file, genbank_file)

    compounds_metanetx, reactions_metanetx, genes_metanetx = map_ec_to_metanetx(metanetx_file, genbank_file)

    draft_genes = set(list(genes_bigg.keys()) + list(genes_modelseed.keys()))

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    with open(output_folder + '/draft.tsv', 'w') as output_file:
        csvwriter = csv.writer(output_file, delimiter='\t')
        csvwriter.writerow(['gene', 'bigg_reactions', 'modelseed_reactions', 'metanetx_reactions'])
        for gene in draft_genes:
            if gene in genes_bigg:
                bigg_reactions = ','.join(genes_bigg[gene])
            else:
                bigg_reactions = ''
            if gene in genes_modelseed:
                modelseed_reactions = ','.join(genes_modelseed[gene])
            else:
                modelseed_reactions = ''
            if gene in genes_metanetx:
                metanetx_reactions = ','.join(genes_metanetx[gene])
            else:
                metanetx_reactions = ''
            csvwriter.writerow([gene, bigg_reactions, modelseed_reactions, metanetx_reactions])

    create_sbml(reactions_bigg, compounds_bigg, output_folder + '/bigg.sbml')
    create_sbml(reactions_modelseed, compounds_modelseed, output_folder + '/modelseed.sbml')
    create_sbml(reactions_metanetx, compounds_metanetx, output_folder + '/metanetx.sbml')