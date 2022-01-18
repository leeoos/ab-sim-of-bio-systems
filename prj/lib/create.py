# This scrip create a LAMMPS input file 
# based on the SBM model given in input

#       --- CLASSES ---
import os
from functools import reduce
import math

# get info of the species
class SpeciesClass:

    def __init__(self, model) -> None:
        self.num_of_species = model.getNumSpecies() 
        self.species = self.__get_species(model)
        self.amounts = self.__compute_amounts(model)   
        self.total_amount = int(reduce(lambda x,y: x+y, self.amounts))
        self.dictionary = self.__make_dict_of_species(model)
        pass

    def __get_species(self, model):
        species = []
        for i in range(self.num_of_species):
            species.append(model.getSpecies(i).getId())
        return species

    def __compute_amounts(self, model):
        amounts = []; scale = True
        for i in range(self.num_of_species):
            specie = model.getSpecies(i)
            value = specie.getInitialConcentration()
            if (math.isnan(value)) : scale = False; break
            else : amounts.append(value)
        if scale:
            # compute a reduced number of moles
            my_n_moli = (reduce(lambda x,y: x+y, amounts))/self.num_of_species * 10
            for i in range(len(amounts)):
                amounts[i] = math.ceil(amounts[i] * my_n_moli)
        else:
            amounts = []
            for i in range(self.num_of_species):
                specie = model.getSpecies(i)
                value = specie.getInitialAmount()
                print(value)
                if (math.isnan(value)) : 
                    print('Error: no initial ammount/concentration given')
                    os._exit(2)
                else :
                    amounts.append(value)
        return amounts

    # entry structure: 'species_id': (atom_id, 'compartment_id', )
    def __make_dict_of_species(self, model):
        dic_of_species = {} ; s_atom_id = 1
        for i in range(self.num_of_species):
            specie = model.getSpecies(i)
            dic_of_species[specie.getId()]= (  s_atom_id,
                                                specie.getCompartment(),
                                                int(self.amounts[i])
                                            ) 
            s_atom_id += 1
        return dic_of_species

# get info of the reactions
class ReactionClass:

    def __init__(self, model) -> None:
        self.model = model
        self.num_of_reactions = model.getNumReactions()
        self.reactions = self.__get_reactions(model)

        __reac_tuple = self.__group_reactants(model)
        self.groups_of_reactants = __reac_tuple[0]
        self.reactants = __reac_tuple[1]

        __prod_tuple = self.__group_products(model)
        self.groups_of_products = __prod_tuple[0]
        self.products = __prod_tuple[1]

        self.combinations = self.__arrange_reactions()
        pass

    def __get_reactions(self, model):
        reactions = []
        for rex in range(self.num_of_reactions):
            reactions.append(model.getReaction(rex).getId())
        return reactions

    def __group_reactants(self, model) :
        reactants = []
        groups_of_reactants = []
        for i in range(self.num_of_reactions):
            rxn = model.getReaction(i)
            rri = []
            for j in range(rxn.getNumReactants()):
                reactant = rxn.getReactant(j).getSpecies()
                rri.append(reactant) ; reactants.append(reactant)
            groups_of_reactants.append(rri)
        return (groups_of_reactants, reactants)

    def __group_products(self, model) :
        products = []
        groups_of_products= []
        for i in range(self.num_of_reactions):
            rxn = model.getReaction(i)
            pri = []
            for j in range(rxn.getNumProducts()):
                reactant = rxn.getProduct(j).getSpecies()
                pri.append(reactant) ; products.append(reactant)
            groups_of_products.append(pri)
        return (groups_of_products, products)

    # complexity = O(m (n!/k!(n-k)!)) 
    # where m = numb of reactions;
    #       n = max numb of reactant per reaction
    #       k = 2
    def __arrange_reactions(self):
        combo = []
        for r_id in range(0, self.num_of_reactions):
            rri = len(self.groups_of_reactants[r_id])
            for i in range(0, rri):
                for j in range(i+1, rri): 
                    combo.append((self.groups_of_reactants[r_id][i],
                                    self.groups_of_reactants[r_id][j]))
        return combo


#           --- FUNCTIONS ---
import os.path
import time
from libsbml import *

#       -- Check SBML file  --
def read_sbml(filename):

    print("\nReading SBML file...")
    time.sleep(1)
    document = readSBML(filename)

    # check for errors in doc
    if document.getNumErrors() > 0:
        print("Encountered the following SBML errors:" )
        document.printErrors()
        os._exit(1)

    # get doc version 
    level = document.getLevel()
    version = document.getVersion()
    short_filename = filename[filename.rfind('/')+1:]

    print("\n"  + "File: " + short_filename
                + " (Level " + str(level) + ", version " + str(version) + ")" )
    time.sleep(1)

    # make a libSBML model to check
    model = document.getModel()

    # check if model is present
    if model is None:
        print("No model present." )
        os._exit(1)

    # check if model is valid
    if model == 1:
        print("No model present." )
        os._exit(1)

    # check if the model is a 
    # Systems Biology Ontology type
    if model.isSetSBOTerm():
        print("      model sboTerm: " + model.getSBOTerm() )
        time.sleep(1)

    return model

#       -- Create a LAMMPS input file --
def make_lmp(**kwargs):

    # optional arguments
    r_seed = kwargs.get('r_seed', 5783)                 # random seed used inside the simulation
    lmp_file_path = kwargs.get('lmp_file_path', None)   # path to dir for lammps input
    sbml_model_file = kwargs.get('sbml_filename', None)   # file name of sbml input model 

    # optional: config some parameters just in case this script is runned indipendently
    if(lmp_file_path == None) :  lmp_file_path = 'in.lmp'

    if (sbml_model_file == None) :
        test = 'test.xml'
        if(os.path.isfile('/home/leeoos/Projects/Tesi/AB-Sim-Of-Bio-Systems/models/'+test)) : 
            sbml_model_file = '/home/leeoos/Projects/Tesi/AB-Sim-Of-Bio-Systems/models/'+test #Alharbi2020
        else: 
            sbml_model_file = input("Insert the path to the a SBML file: ")

    # read sbml input doc and check for errors 
    # if errors then return exit code 1
    model = read_sbml(sbml_model_file)
    short_filename = sbml_model_file[sbml_model_file.rfind('/')+1:-4]

    print("\nCreating lammps file as in.lmp ...")
    time.sleep(1)
    
    # new objects of model's classes such as SpeciesClass, ReactionClass ...
    S = SpeciesClass(model)
    R = ReactionClass(model)

    print()
    print(S.dictionary)
    print(R.combinations)

    for i in range(R.num_of_reactions):
        print(R.reactions[i], ": ", '\n')
        print("         ", end=' ')
        for react in R.groups_of_reactants[i] :
            print(react,"("+str(S.dictionary[react][0])+")",  end=' ')
        print("     -->       ", end='')
        for prod in R.groups_of_products[i] :
            print(prod,"("+str(S.dictionary[prod][0])+")",  end=' ')
        print()
        print("\n")
    

    with open(lmp_file_path, 'w') as f:
        f.write('# Agent Based Simulation Of Biological Systems\n\n')
    
        # SET UP OF INPUT VARIABLES
        set_up =[
        "\n#       --- SET UP OF INPUT VARIABLES ---\n",

        "# duration : int = N number of steps for the current run",
        "variable time_value index 50  # default value: 50 loop * 100 steps",
        "variable loop_len equal ${time_value}\n",

        "# atoms : int = N number of atoms of each type to generate ",
        "# uncomment this lines and substitue the number of atoms",
        "# in the create_atoms command with this variable ",
        "# to modify the initial amount of atoms in the simulation;",
        "# prj/run.sh -h for more info",
        "# Note: doing so the initial amount will be the same for all atoms types",
        "#variable num_atoms index 5",
        "#variable atoms equal ${num_atoms}\n"]

        f.writelines(["%s\n" % item  for item in set_up])

        if (S.total_amount == 0 ) : expeted_types = S.num_of_species + 1
        else : expeted_types = S.num_of_species
    
        # SIMULATION BOX PROPERTIES
        sim_box = [
        "\n#       --- SIMULATION BOX PROPERTIES ---\n",

        "# sym measure units and atoms style",
        "units       lj",
        "atom_style  full \n",

        "# box dimension, boudaries and structure",
        "dimension   3  ",
        "boundary    f f f ",
        "lattice     fcc 3.52\n",

        "# flag2 = on or off for bonded interactions",
        "newton on off\n",

        "# define simulation box",
        "region      box block 0 30 0 30 0 30",
        "create_box  "+str(expeted_types)+" box  bond/types " + str(len(R.reactants)) + " extra/bond/per/atom 100\n",

        "# create simulation walls",
        "fix xwalls all wall/reflect xlo EDGE xhi EDGE",
        "fix ywalls all wall/reflect ylo EDGE yhi EDGE",
        "fix zwalls all wall/reflect zlo EDGE zhi EDGE\n"]

        f.writelines(["%s\n" % item  for item in sim_box])

        # AGENTS PROPRETIES AND FORCE FIELDS
        f.write("\n#       --- AGENTS PROPRETIES AND FORCE FIELDS ---\n")
        f.write("# creation of atoms of types in randoms spots inside the box\n")
        
        for k in list(S.dictionary.keys()):
            types = ""
            if k in R.reactants : types = types + "     | reactant |"
            if k in R.products : types = types + "      | product  |"
            f.write("create_atoms"+"    "+ str(S.dictionary[k][0]) +" random "+ str(S.dictionary[k][2]) +
                        " "+ str(r_seed) + " box  # "+ k + types)
                        
            r_seed = r_seed + 1
            f.write("\n")

        if (S.total_amount == 0): 
            f.write("create_atoms"+"    "+ str(S.num_of_species+1) +" random 5 "
                        + str(r_seed) + " box  # ghost atom to inizialize velocity")
            r_seed = r_seed + 1
            f.write("\n")

        f.write("\n# atoms mass\n")
        for i in list(S.dictionary.values()) :
            f.write("mass " + str(i[0]) + " 1.0")
            f.write("\n")
        if (S.total_amount == 0 ) : 
            f.write("mass " + str(S.num_of_species+1) + " 1.0")
            f.write("\n")

        f.write("\n# assing atoms to agents groups\n")
        
        for i in range (0, R.num_of_reactions):
            types = ""
            new_group = "group agents" + str(i+1) + " type "
            if (R.groups_of_reactants[i] == []) : pass
            else :
                for k in R.groups_of_reactants[i]:
                    types = types + str(S.dictionary[k][0]) + " "
                f.write(new_group + types)
                f.write("\n")

        if (S.total_amount == 0) :
            f.write("group to_dump empty\n")
 
        ag_prop = [
        "\n# force fields style and coefficient",
        "pair_style zero 5.0",
        "pair_coeff * *\n",
        
        "# bond style and coefficients",
        "bond_style  harmonic",
        "bond_coeff * 100 1.1\n"]

        f.writelines(["%s\n" % item  for item in ag_prop])
        
        # SIMULATION 
        sim1 = [
        " \n#      --- SIMULATION ---\n",

        "# set time steps ",
        "timestep 0.001   # seconds\n",
        
        "# This command sets parameters that affect",
        "# the building of pairwise neighbor lists",
        "neighbor 0.001 bin",
        "neigh_modify every 10 delay 100\n",
        
        "# print thermodinamic inf every N timesteps",
        "thermo 100 \n",
        
        "# fix ID group-ID bond/create Nevery itype jtype Rmin bondtype keyword values",
        "# this fix will attempt to create new bond btw atoms of type i and j every N timestep"]

        f.writelines(["%s\n" % item  for item in sim1])

        for i, j in zip(R.combinations, range(len(R.combinations))):
            f.write("fix bond"+str(j+1)+" all bond/create 10 "+ str(S.dictionary[i[0]][0]) +" "+ 
                                str(S.dictionary[i[1]][0])+ " 1.0 "+ str(j+1) +" prob 0.5 " + str(r_seed))
            r_seed = r_seed + 1
            f.write("\n")
        
        sim2 = [
        "\n# set velocity for all atoms",
        "velocity all create 300.0 4928459 rot yes dist gaussian \n",
        
        "# perform plain time integration ",
        "# to update position and velocity",
        "# and simulate Brownioan motion",
        "fix 1 all nve",
        "fix 2 all langevin 300.0 300.0 10.0 904297\n",
        
        "# compute if atoms has a bonds",
        "# and total number of bonds btw all atoms"]

        f.writelines(["%s\n" % item  for item in sim2])

        thermo_style1 = "thermo_style custom step temp pe "
        thermo_style2 = "thermo_style custom step temp pe "
        counter = "variable counter"

        if (S.total_amount == 0) :
            f.write("compute t   all property/atom type\n")
            to_dump = "to_dump"
            righttype = ("# righttype : boolean = true if atom I is of type between 1 and " 
                        + str(S.num_of_species) + 
                        "\nvariable righttype atom 'c_t < "+ str(S.num_of_species+1) +"'\n")
        else: 
            righttype = ""
            to_dump = "all"

        f.write("compute hb0 all property/atom nbonds \n\n")

        for i in range(1, R.num_of_reactions+1):
            if (R.groups_of_reactants[i-1] == []) : pass
            else:
                f.write("compute hb"+str(i)+" agents"+str(i)+" property/atom nbonds\n")
                f.write("compute cb"+str(i)+" agents"+str(i)+" reduce sum c_hb"+str(i)+"\n")

                thermo_style1 = thermo_style1 + "c_cb"+str(i) +" "
                counter = counter + str(i) + " equal ceil(c_cb" + str(i) + ")\nvariable counter"
                thermo_style2 = thermo_style2 + "v_counter" + str(i) + " "

                f.write("\n")

        sim3 = [
        "\n# this lines are necessary to insure that the “hasbond” and 'newatoms' ",
        "# variables are current when the group command invokes it.",          
        thermo_style1,
        "run 0",

        "\n# hasbond : boolean = true if atom I has a bond with atom J",
        "variable hasbond atom 'c_hb0 > 0.0'",

        righttype,

        "\n# counter : int = N total number of bonds in the sim",
        counter[:counter.rfind('\n')],
        
        "\n# print themo info every timestep",
        thermo_style2,
        
        "\n# dumps atoms information",
        "dump 1 "+ to_dump + " custom 10 dump."+short_filename+" id x y z type \n"]

        f.writelines(["%s\n" % item  for item in sim3])

        # LOOP 
        loop1 = [
        "\n#       --- LOOP---\n",
        
        "label loop",
        "variable step loop ${loop_len}   # loop length\n",
        
        "# create new atoms only if new bonds have been made", 
        "# the num of new atoms is linked to the number of new bonds as follow:",
        "# fix ID group-ID deposit N type M seed keyword values"]

        f.writelines(["%s\n" % item  for item in loop1])

        for i in range(len(R.groups_of_products)):
            if (R.groups_of_products[i] == []) : pass
            else :
                for p in R.groups_of_products[i] :
                    if (R.groups_of_reactants[i] == []):
                        f.write("fix deposit"+str(i+1)+" all deposit 1 "+str(S.dictionary[p][0])
                                    +" 1 5748 region box near 2.0\n")
                    else:
                        f.write("variable newatoms"+str(i+1)+" equal floor(${counter"+str(i+1)+"}/2)\n")
                        f.write("if '${counter"+str(i+1)+"} > 0' then &\n")
                        f.write("'fix deposit"+str(i+1)+" all deposit ${newatoms"+str(i+1)+"} "+str(S.dictionary[p][0])
                                    +" 1 5748 region box near 2.0' \n")
                        f.write("\n")

        if(S.total_amount == 0):
            to_dump = ("# assing all atoms of the right kind to the dump group\n"
                    +"group to_dump dynamic all every 1 var righttype \n")
        else: to_dump = ""
        
        loop2 = [
        "\n# assing all atoms that have a bond to the garbage group",
        "group garbage dynamic all every 1 var hasbond",

        to_dump,
        
        "# append new values on dump file",
        "dump_modify 1 append yes\n",
        
        "# perform n steps in loop",
        "run 100\n",
    
        "# delate all atoms in garbage",
        "delete_atoms group garbage bond yes mol yes compress no\n",
        
        "# jump to loop lable until step > 0 ",
        "next step",
        "jump SELF loop\n",
        
        "# end of loop",
        "label break\n",
        
        "# check on input variables",
        "variable duration equal ${loop_len}*100",
        "print ''",
        "print 'Starting Atoms: "+ str(S.total_amount) +" ' ",
        "print 'Duration: ${duration}'",
        "print 'ALL DONE' \n"]

        f.writelines(["%s\n" % item  for item in loop2])

        f.flush()
        os.fsync(f)

if __name__ == '__main__':
    make_lmp()