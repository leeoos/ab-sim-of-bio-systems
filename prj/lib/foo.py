 
#from libsbml import *

import libsbml



def make_lmp(**kwargs):

    sbml_test_file = '/home/leeoos/Projects/Tesi/AB-Sim-Of-Bio-Systems/resources/sbmlex/test.xml'

    reader = libsbml.SBMLReader()
    document = reader.readSBML(sbml_test_file)
    model = document.getModel()

    groups_of_reactants = []
    for r in range(model.getNumReactions()):
        rxn = model.getReaction(r)
        print(rxn)
        nor = rxn.getNumReactants()
        rri = []
        for m in range(nor):
            reactant = rxn.getReactant(m).getSpecies()
            print(r, m, reactant)
            rri.append(reactant)
        groups_of_reactants.append(rri)

    species = []
    for r in range(model.getNumSpecies()):
        rxn = model.getSpecies(r)
        species.append(rxn.getId())
 
    print(groups_of_reactants)
    print(species)
 

if __name__ == '__main__':
    make_lmp()