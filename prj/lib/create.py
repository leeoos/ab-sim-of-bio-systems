# This scrip create a LAMMPS input file 
# based on the SBM model given in input

class SpeciesClass:

    def __init__(self, model) -> None:
        self.model = model.getNumSpecies()
        self.num_of_species = model.getNumReactions()
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
        self.groups_of_reactants = self.group_reactants()
        self.groups_of_products = self.group_products()
        self.reactants = ['s1', 's3', 's4', 's2', 's5']
        self.products = ['s6', 's7']
        self.combinations = self.arrange_reactions()
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

    # make_lmp optional arguments
    r_seed = kwargs.get('r_seed', 5783)                 # random seed used inside the simulation
    lmp_file_path = kwargs.get('lmp_file_path', None)   # path to dir for lammps input
    sbml_filename = kwargs.get('sbml_filename', None)   # file name of sbml input model 

    # optional: config some parameters just in case this script is runned indipendently
    if(lmp_file_path == None) :  lmp_file_path = 'in.lmp'

    if (sbml_filename == None) :
        if(os.path.isfile('/home/leeoos/Projects/Tesi/AB-Sim-Of-Bio-Systems/resources/sbmlex/test.xml')) : 
            sbml_filename = '/home/leeoos/Projects/Tesi/AB-Sim-Of-Bio-Systems/resources/sbmlex/test.xml'
        else: 
            sbml_filename = input("Insert the path to the a SBML file: ")

    # read sbml input doc and check for errors 
    # if errors then return exit code 1
    model = read_sbml(sbml_filename)

    if model == 1:
        os._exit(1)

    print("\nCreating lammps file as in.lmp ...")
    time.sleep(1)
    
    S = SpeciesClass(model)
    R = ReactionClass(model)

    with open(lmp_file_path, 'w') as f:
        f.write('# Agent Based Simulation Of Biological Systems\n\n')
    
        # SET UP OF INPUT VARIABLES
        set_up =[
        "\n#       --- SET UP OF INPUT VARIABLES ---\n",

        "# rnseed : int = seed for random numbers",
        "variable rnseed index 10",
        "variable probability equal random(0,1,${rnseed})\n",

        "# duration : int = N number of steps for the current run",
        "variable time_value index 50  #default value",
        "variable duration equal ${time_value}\n",

        "# atoms : int = N number of atoms of each type to generate ",
        "variable num_atoms index 5",
        "variable atoms equal ${num_atoms}\n"]

        f.writelines(["%s\n" % item  for item in set_up])
    
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
        "create_box  "+str(S.num_of_species)+" box  bond/types " + str(len(R.reactants)) + " extra/bond/per/atom 100\n",

        "# create simulation walls",
        "fix xwalls all wall/reflect xlo EDGE xhi EDGE",
        "fix ywalls all wall/reflect ylo EDGE yhi EDGE",
        "fix zwalls all wall/reflect zlo EDGE zhi EDGE\n"]

        f.writelines(["%s\n" % item  for item in sim_box])

        # AGENTS PROPRETIES AND FORCE FIELDS
        f.write("\n#       --- AGENTS PROPRETIES AND FORCE FIELDS ---\n")
        f.write("# creation of atoms of types in randoms spots inside the box\n")
        
        for k in list(S.dictionary.keys()):
            if k in R.products:
                f.write("create_atoms"+"    "+ str(S.dictionary[k][0]) +" random 0 "+ str(r_seed) + " box  # "+ k 
                        +": product ")
            else:
                f.write("create_atoms"+"    "+ str(S.dictionary[k][0]) +" random ${atoms} "+ str(r_seed) + " box  # "+k
                        +": reactant ")
                        
            r_seed = r_seed + 1
            f.write("\n")

        f.write("\n# atoms mass\n")
        for i in list(S.dictionary.values()) :
            f.write("mass " + str(i[0]) + " 10.948")
            f.write("\n")

        f.write("\n# assing atoms to agents groups\n")
        
        for i in range (0, R.num_of_reactions):
            types = ""
            new_group = "group agents" + str(i+1) + " type "
            for k in R.groups_of_reactants[i]:
                types = types + str(S.dictionary[k][0]) + " "
            f.write(new_group + types)
            f.write("\n")

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
        "# this fix will attempt to create new bond btw atoms of ",
        "# type 1 and 2 every Nevery timestep"]

        f.writelines(["%s\n" % item  for item in sim1])

        for i, j in zip(R.combinations, range(0, len(R.combinations))) :
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

        f.write("compute hb0 all property/atom nbonds \n\n")

        for i in range(1, R.num_of_reactions+1):
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

        "\n# counter : int = N total number of bonds in the sim",
        counter[:counter.rfind('\n')],
        
        "\n# print themo info every timestep",
        thermo_style2,
        
        "\n# dumps atoms information",
        "dump 1 all custom 10 dump.out id x y z type \n"]

        f.writelines(["%s\n" % item  for item in sim3])

        # LOOP 
        loop1 = [
        "\n#       --- LOOP---\n",
        
        "label loop",
        "variable step loop ${duration}   # loop length\n",
        
        "# create new atoms only if new bonds have been made", 
        "# the num of new atoms is linked to the number of new bonds as follow:",
        "# fix ID group-ID deposit N type M seed keyword values"]

        f.writelines(["%s\n" % item  for item in loop1])

        for p, i in zip(R.products, range(1, R.num_of_reactions +1)):
            f.write("variable newatoms"+str(i)+" equal floor(${counter"+str(i)+"}/2)\n")
            f.write("if '${counter"+str(i)+"} > 0' then &\n")
            f.write("'fix depositatoms all deposit ${newatoms"+str(i)+"} "+str(S.dictionary[p][0])+" 1 5748 region box near 2.0' \n")
            f.write("\n")

        "'fix depositatoms all deposit ${newatoms} 6 1 5748 region box near 2.0' ",
        "\n",
        
        loop2 = [
        "\n# assing all atoms that have a bond to the garbage group",
        "group garbage dynamic all every 1 var hasbond\n",
        
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
        "variable total_atoms equal ${num_atoms}*2",
        "print ''",
        "print 'Starting Atoms: ${total_atoms}' ",
        "print 'Duration: ${duration}'",
        "print 'ALL DONE' \n"]

        f.writelines(["%s\n" % item  for item in loop2])

        f.flush()
        os.fsync(f)

if __name__ == '__main__':
    make_lmp()