# /bin/python3
# This scrip create a LAMMPS input file 
# based on the SBM model given in input

#       --- CLASSES ---
from functools import reduce
from decimal import Decimal
from random import randint
import math
import os

# get info about the
class CompartmentClass:

    def __init__(self,model) -> None:
        if (model.getNumCompartments()) > 1 :
            print(  "\nWarning: all the agents will be in the same compartment. "
                    +"If agents of the same species are in different compartment "
                    +"they will be considered as two diffrent kind of agents "
                    +"this means that the system may behave unpredictably\n"
                )
        self.csize = int(model.getCompartment(0).getSize())
        pass

# get info about the species
class SpeciesClass:

    def __init__(self, model) -> None:
        self.num_of_species = model.getNumSpecies() 
        self.initial_atoms = self.__compute_initial_atoms(model)   
        self.total_initial_atoms = int(reduce(lambda x,y: x+y, self.initial_atoms))
        self.dictionary = self.__make_dict_of_species(model)
        pass

    def __fexp(self, number):
        (sign, digits, exponent) = Decimal(number).as_tuple()
        return len(digits) + exponent - 1

    def __compute_initial_atoms(self, model):
        C = CompartmentClass(model)
        initial_atoms = []; scale = True
        for i in range(self.num_of_species):
            specie = model.getSpecies(i)
            value = specie.getInitialConcentration()
            if (math.isnan(value)) : scale = False; break
            else : initial_atoms.append(math.ceil(value* pow(10,-(self.__fexp(value)))) * C.csize)
        if not(scale):
            initial_atoms = []
            for i in range(self.num_of_species):
                specie = model.getSpecies(i)
                value = specie.getInitialAmount()
                if (math.isnan(value)) : 
                    print('Error: no initial ammount/concentration given')
                    os._exit(3)
                else :
                    initial_atoms.append(value)
        return initial_atoms

    # entry structure: 'species_id':[
    #                                   0: atom_id, 
    #                                   1: 'compartment_id',
    #                                   2: initial_amount, 
    #                                   3: groth_rate, 
    #                                   4: deposit_counter 
    #                               ]
    def __make_dict_of_species(self, model):
        dic_of_species = {} ; s_atom_id = 1
        for i in range(self.num_of_species):
            specie = model.getSpecies(i)
            initial_atoms_amount = int(self.initial_atoms[i])
            growth_rate = randint(1, 10)
            deposit_counter = 0
            dic_of_species[specie.getId()] = [  
                                                s_atom_id,
                                                specie.getCompartment(),
                                                initial_atoms_amount,
                                                growth_rate,
                                                deposit_counter
                                            ] 
            s_atom_id += 1
        return dic_of_species

# get info about the reactions
class ReactionClass:

    def __init__(self, model) -> None:
        self.reactions = self.__get_reaction_info(model)
        self.num_of_reactions = len(self.reactions)
        self.combinations = self.__arrange_reactions()
        pass

    def __get_reaction_info(self, model):
        self.reactants = set(); self.products = set()
        reactions = {}
        for i in range(model.getNumReactions()):
            rxn = model.getReaction(i)
            rri = [] # reatants of reaction i
            for j in range(rxn.getNumReactants()):
                reactant = rxn.getReactant(j).getSpecies()
                rri.append(reactant) ; self.reactants.add(reactant)
            pri = [] # products of reaction i
            for k in range(rxn.getNumProducts()):
                product = rxn.getProduct(k).getSpecies()
                pri.append(product) ; self.products.add(product)
            reactions[rxn.getId()] = (rri, pri)
        return reactions

    # complexity = O(m (n!/k!(n-k)!)) 
    # where m = numb of reactions;
    #       n = max numb of reactant per reaction
    #       k = 2
    def __arrange_reactions(self):
        combo = []
        for r_id in list(self.reactions.keys()): 
            rri = len(self.reactions[r_id][0])
            for i in range(0, rri):
                for j in range(i+1, rri): 
                    combo.append((self.reactions[r_id][0][i],
                                    self.reactions[r_id][0][j]))
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

    # analyze the sbml document and dump the info 
    os.system('rm sbml.analysis 2> /dev/null')
    analysis = 'sbml.analysis'
    with open(analysis, 'w') as m:
        m.write("Analysis of SBML file: "+ short_filename +"\n")
        m.write("\n\nSpecies \n\n")
        m.write("{:<20s} {:<20s} {:<25s} {:<30s}\n\n".format("Species id", "Atom id", "Compartment", "Amount")) 
        for key, value in S.dictionary.items():
            specie = str(key); atomid = str(value[0]); comp = str(value[1]); amo = str(value[2])
            m.write("{:<20s} {:<20s} {:<25s} {:<30s}\n".format(specie, atomid, comp, amo))
            m.write("\n")
        m.write("\nReactions Map \n\n")
        react = [(R.reactions[r_id][0], R.reactions[r_id][1]) for r_id in list(R.reactions.keys())]
        for r_id, i in zip(list(R.reactions.keys()), range(R.num_of_reactions)):
            m.write("-"+ r_id +": \n\n")
            m.write("\t\t{:<20s}  ->  {:<20s}\n".format(str(react[i][0]), str(react[i][1])))
            m.write("\n\n")
        m.write("\n\n")
        m.flush()
        os.fsync(m)
        
    os.system('rm '+ lmp_file_path +' 2> /dev/null')
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

        if (S.total_initial_atoms == 0 ) : expeted_types = S.num_of_species + 1
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
        "create_box  "+str(expeted_types)+" box  bond/types 1 extra/bond/per/atom 100\n", # ATTENZIONE ALL'EXTRA BOND TYPE

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

        if (S.total_initial_atoms == 0): 
            f.write("create_atoms"+"    "+ str(S.num_of_species+1) +" random 5 "
                        + str(r_seed) + " box  # ghost atom to inizialize velocity")
            r_seed = r_seed + 1
            f.write("\n")

        f.write("\n# atoms mass\n")
        for i in list(S.dictionary.values()) :
            f.write("mass " + str(i[0]) + " 1.0")
            f.write("\n")
        if (S.total_initial_atoms == 0 ) : 
            f.write("mass " + str(S.num_of_species+1) + " 1.0")
            f.write("\n")

        f.write("\n# assing atoms to agents groups\n")
        mortals = set() ; agents_counter = 1
        for r_id in list(R.reactions.keys()): #range (R.num_of_reactions):
            new_group = "group agents" + str(agents_counter) + " type " ; types = ""
            only_one = (len(R.reactions[r_id][0]) ==  1)
            if (R.reactions[r_id][0] == []) : pass
            else :
                for i in R.reactions[r_id][0]:
                    if (only_one ) : mortals.add(str(S.dictionary[i][0]))
                    else: types = types + str(S.dictionary[i][0]) + " "
                if (not only_one) : 
                    f.write(new_group + types)
                    f.write("\n")
            agents_counter += 1

        if (mortals != set()) :
            mortals = ' '.join(mortals)
            f.write("group mortals type "+ mortals +" \n") 

        if (S.total_initial_atoms == 0) :
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
        # EDIT TYPE OF BOND AND PROB 
        for i in R.combinations :
            b1 = str(S.dictionary[i[0]][0]) ; b2 = str(S.dictionary[i[1]][0])
            f.write("fix bond_"+ b1 +"_"+ b2 +" all bond/create 10 "+ 
                        b1 +" "+ b2 +" 1.0 1 prob 0.5 " + str(r_seed))
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

        if (S.total_initial_atoms == 0) :
            f.write("compute t   all property/atom type\n")
            to_dump = "to_dump"
            righttype = ("# righttype : boolean = true if atom I is of type between 1 and " 
                        + str(S.num_of_species) + 
                        "\nvariable righttype atom 'c_t < "+ str(S.num_of_species+1) +"'\n")
        else: 
            righttype = ""
            to_dump = "all"

        f.write("compute hb0 all property/atom nbonds \n\n")

        cc = 1  #compute counter
        for r_id in list(R.reactions.keys()): #range(R.num_of_reactions):
            if (len(R.reactions[r_id][0]) <= 1) :  pass
            else:
                f.write("compute hb"+str(cc)+" agents"+str(cc)+" property/atom nbonds\n")
                f.write("compute cb"+str(cc)+" agents"+str(cc)+" reduce sum c_hb"+str(cc)+"\n")

                thermo_style1 = thermo_style1 + "c_cb"+str(cc) +" "
                counter = counter + str(cc) + " equal ceil(c_cb" + str(cc) + ")\nvariable counter" 
                thermo_style2 = thermo_style2 + "v_counter" + str(cc) + " "

                f.write("\n")
            cc += 1

        if (counter == "variable counter" ) : counter = "# no bonds detected "

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
        "dump 1 "+ to_dump + " custom 10 dump."+ short_filename +".out id x y z type \n"]

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

        cc = 1  #compute counter
        for r_id in list(R.reactions.keys()):
            deposit = []
            if (R.reactions[r_id][1] == []) :  pass
            else :
                for p in R.reactions[r_id][1] :
                    atom = str(S.dictionary[p][0]) +"_"+ str(S.dictionary[p][4])
                    S.dictionary[p][4] += 1
                    p_id = str(S.dictionary[p][0])
                    growth = S.dictionary[p][3]
                    if (len(R.reactions[r_id][0]) <= 1) :
                        f.write("fix deposit"+ atom +" all deposit 1 "+ p_id +" "
                                    + str(growth) +" 5748 region box near 2.0\n")
                        deposit.append(0)
                    else:
                        deposit.append("'fix deposit"+ atom +" all deposit ${newatoms"+ str(cc) 
                                        +"} "+ p_id +" "+ str(growth) +" 5748 region box near 2.0' &")

                if (deposit[0] == 0) :
                    f.write("\n") ; pass
                else:
                    f.write("variable newatoms"+str(cc)+" equal floor(${counter"+str(cc)+"}/2)\n")
                    f.write("if '${counter"+str(cc)+"} > 0' then &\n")
                    deposit[len(deposit)-1] = deposit[len(deposit)-1][:-1]
                    f.writelines(["%s\n" % item  for item in deposit])
                    f.write("\n")
            cc += 1

        if(S.total_initial_atoms == 0):
            to_dump = ("# assing all atoms of the right kind to the dump group\n"
                    +"group to_dump dynamic all every 1 var righttype \n")
        else: to_dump = ""


        if (mortals != set()) :
            inhibition = str(randint(1,10))
            mortals = "# delate M atoms in mortals every N timestamps \n"
            mortals = mortals + "# fix ID group-ID evaporate N M region-ID seed \n"
            mortals = mortals + "fix death mortals evaporate "+ inhibition +" 1 box "+ str(r_seed) +"\n"
                            
        
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

        mortals,
        
        "# jump to loop lable until step > 0 ",
        "next step",
        "jump SELF loop\n",
        
        "# end of loop",
        "label break\n",
        
        "# check on input variables",
        "variable duration equal ${loop_len}*100",
        "print ''",
        "print 'Starting Atoms: "+ str(S.total_initial_atoms) +" ' ",
        "print 'Duration: ${duration}'",
        "print 'ALL DONE' \n"]

        f.writelines(["%s\n" % item  for item in loop2])

        f.flush()
        os.fsync(f)

        print("\nALL DONE\n")

import logging
import sys
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

if __name__ == '__main__':
    make_lmp()