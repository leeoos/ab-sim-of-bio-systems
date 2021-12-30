# This scrip create a LAMMPS input file 
# based on the SBM model given in input

class SpeciesClass:

    def __init__(self, model) -> None:
        self.num_of_species = model.getNumSpecies() 
        self.dictionary = self.make_dict()
        pass

    # DA RIFINIRE
    def make_dict(self):
        reactants = ['s1', 's3', 's4', 's2', 's5']
        products = ['s6', 's7']
        r_id_dic = {} ; r_value1 = 1
        for i in reactants:
            r_id_dic[i] = (r_value1,0) ; r_value1 += 1
        for i in products:
            r_id_dic[i] = (r_value1,0) ; r_value1 += 1
        return r_id_dic


class ReactionClass:

    def __init__(self, model) -> None:
        self.model = model
        self.num_of_reactions = model.getNumReactions()
        self.list_of_reactions = model.getListOfReactions()
        self.groups_of_reactants = self.group_reactants()
        self.groups_of_products = self.group_products()
        self.reactants = ['s1', 's3', 's4', 's2', 's5']
        self.products = ['s6', 's7']
        #self.combinations = self.arrange_reactions()
        pass

    def group_reactants(self) :
        rr1 = ['s1', 's4', 's3']
        rr2 = ['s2', 's5']
        return [rr1, rr2]

    def group_products(self) :
        rr1 = ['s6']
        rr2 = ['s7']
        return [rr1, rr2]

    # complexity = O(m (n!/k!(n-k)!)) 
    # where m = numb of reactions;
    #       n = max numb of reactant per reaction
    #       k = 2
    def arrange_reactions(self):
        combo = []
        for r_id in range(0, self.num_of_reactions):
            rri = len(self.groups_of_reactants[r_id])
            for i in range(0, rri):
                for j in range(i+1, rri): 
                    combo.append((self.groups_of_reactants[r_id][i],
                                    self.groups_of_reactants[r_id][j]))
        return combo

#       -- Create LAMMPS input file --
import os
import os.path
import time
import re
from libsbml import *

#       -- Read SBML file  --
def read_sbml(filename):

    print("\nReading SBML file...")
    time.sleep(1)
    document = readSBML(filename)

    #       -- Checks for errors --
    if document.getNumErrors() > 0:
        print("Encountered the following SBML errors:" )
        document.printErrors()
        return 1

    #       -- Get doc version --
    level = document.getLevel()
    version = document.getVersion()
    short_filename = filename[filename.rfind('/')+1:]

    print("\n"  + "File: " + short_filename
                + " (Level " + str(level) + ", version " + str(version) + ")" )
    time.sleep(1)

    model = document.getModel()

    if model is None:
        print("No model present." )
        return 1

    if model.isSetSBOTerm():
        print("      model sboTerm: " + model.getSBOTerm() )
        time.sleep(1)

    return model

#       -- Create a LAMMPS input file --
def make_lmp(**kwargs):

    # optional: config some parameters just in case this script is runned indipendently
    lmp_file_path = 'in.lmp'

    sbml_filename = '/home/leeoos/Projects/Tesi/AB-Sim-Of-Bio-Systems/resources/sbmlex/Alharbi2020.xml'


    # read sbml input doc and check for errors 
    # if errors then return exit code 1
    model = read_sbml(sbml_filename)

    if model == 1:
        os._exit(1)

    print("\nCreating lammps file as in.lmp ...")
    time.sleep(1)
    
    S = SpeciesClass(model)
    R = ReactionClass(model)
    
    print()
    print(S.num_of_species)
    print(R.num_of_reactions)
    print

if __name__ == '__main__':
    make_lmp()