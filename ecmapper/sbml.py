from cobra import Model, Reaction, Metabolite
from cobra.io import write_sbml_model

def create_sbml(reaction_genes, reactions_metabolites, output_file):
    model = Model()
    metabolites_created = []

    for reaction_id in reaction_genes:
        reaction = Reaction(reaction_id)
        reaction.name = reaction_id
        reaction_metabolites = {}
        for reactant_id in reactions_metabolites[reaction_id][0]:
            if reactant_id not in metabolites_created:
                reactant_metabolite = Metabolite(reactant_id, compartment='c')
                reaction_metabolites[reactant_metabolite] = -1.0
        for product_id in reactions_metabolites[reaction_id][1]:
            if product_id not in metabolites_created:
                product_metabolite = Metabolite(product_id, compartment='c')
                reaction_metabolites[product_metabolite] = 1.0
        reaction.add_metabolites(reaction_metabolites)
        reaction.notes['GENE_ASSOCIATION'] = '(' + ' or ' .join(reaction_genes[reaction_id]) + ')'
        model.add_reactions([reaction])

    write_sbml_model(model, output_file)