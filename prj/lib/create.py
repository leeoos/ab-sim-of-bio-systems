# This scrip create a LAMMPS input file 
# based on the SBM model given in input

class ReactionValues:

    def __init__(self, model) -> None:
        self.model = model
        self.num_of_reactions = model.getNumReactions()
        self.reactants = self.get_all_reactant()
        self.products = self.get_all_products()
        self.dictionary = self.make_dict()
        pass

    def get_all_reactant(self) :
        rr1 = ['s1', 's4', 's3']
        rr2 = ['s2', 's5']
        return [rr1, rr2]

    def get_all_products(self):
        pr1 = ['s6']
        pr2 = ['s7']
        return [pr1, pr2]

    def make_dict(self):
        r_dic = {} ; r_value = 1
        for k in range(0,self.num_of_reactions):
            for i in self.reactants[k]:
                r_dic[i] = r_value ; r_value += 1
            for i in self.products[k]:
                r_dic[i] = r_value ; r_value += 1
        return r_dic


    def arrange_reaction(self, r_id) :
        foo = []
        copules = []
        num_of_r = len(self.reactants[r_id])
        iterator = range(1, ((num_of_r * (num_of_r -1)) +1 ))
        for i in iterator :
            if (i % num_of_r == 0) :
                foo.append(num_of_r)
            else :
                foo.append(i % num_of_r)

        for k in foo:
            copules.append(self.dictionary[self.reactants[r_id][k-1]])

        copules = [copules[i:i + 2] for i in range(0, len(copules), 2)]
        return copules

#       -- Create LAMMPS input file --

import os
import os.path
import time
from libsbml import *

#       -- Read SBML file  --
def read_sbml(filename):

    #print("\nSelect a SBLM Model")
    #filename = get_file()

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

    r_seed = kwargs.get('r_seed', 5783)
    lmp_file_path = kwargs.get('lmp_file_path', None)
    sbml_filename = kwargs.get('sbml_filename', None)


    if (lmp_file_path == None) : 
        os.system('mkdir -p ./lammps/')
        lmp_file_path = 'lammps/in.lmp'

    if(sbml_filename == None) : sbml_filename = '/home/leeoos/Projects/Tesi/AB-Sim-Of-Bio-Systems/sbmlex/test.xml'
        #sbml_filename = input("Insert the path to the a SBML file: ")


    if (lmp_file_path == None) : 
        os.system('mkdir -p ./lammps/')
        lmp_file_path = 'lammps/in.lmp'

    if(sbml_filename == None) : sbml_filename = input("Insert the path to the a SBML file: ")


    model = read_sbml(sbml_filename)

    if model == 1:
        os._exit(1)

    print("\nCreating lammps file as in.lmp ...")
    time.sleep(1)

    # info
    species = model.getNumSpecies()

    r = ReactionValues(model)
    reactants = [item for sublist in r.reactants for item in sublist]
    print(reactants)
    products = [item for sublist in r.products for item in sublist]
    print(products)
    r_dictionary = r.dictionary
    print(r.dictionary)
    print(r.arrange_reaction(0))
    print(r.arrange_reaction(1))  


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
        "create_box  "+str(species)+" box  bond/types " + str(len(r.reactants)) + " extra/bond/per/atom 100\n",

        "# create simulation walls",
        "fix xwalls all wall/reflect xlo EDGE xhi EDGE",
        "fix ywalls all wall/reflect ylo EDGE yhi EDGE",
        "fix zwalls all wall/reflect zlo EDGE zhi EDGE\n"]

        f.writelines(["%s\n" % item  for item in sim_box])

        # AGENTS PROPRETIES AND FORCE FIELDS
        f.write("\n#       --- AGENTS PROPRETIES AND FORCE FIELDS ---\n")
        f.write("# creation of atoms of types in randoms spots inside the box\n")
        
        for k in list(r_dictionary.keys()):
            if k in products:
                f.write("create_atoms" + "    " + str(r_dictionary[k]) + " random 0 " + str(r_seed) + " box  #product")
            else:
                f.write("create_atoms" + "    " + str(r_dictionary[k]) + " random ${atoms} " + str(r_seed) + " box")
            r_seed = r_seed + 1
            f.write("\n")

        f.write("\n# atoms mass\n")
        for i in list(r_dictionary.values()) :
            f.write("mass " + str(i) + " 10.948")
            f.write("\n")

        f.write("\n# assing atoms to cerian groups\n")
        for k in list(r_dictionary.keys()) :
            if k in reactants:
                f.write("group g" + str(r_dictionary[k]) + " type "+str(r_dictionary[k]))
            f.write("\n")
        j=1
        for i in r.reactants:
            group_name = "group agents"+str(j)+" union "
            j += 1
            for k in i: group_name = group_name + "g"+str(r_dictionary[k]) + " "
            f.write(group_name)
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

        j = 1
        for i in range(0, r.num_of_reactions):
            for k in r.arrange_reaction(i) :
                f.write("fix bond"+str(j)+" all bond/create 10 "+ str(k[0]) +" "+ str(k[1])+ " 1.0 "+ str(i+1) +" prob 0.5 " + str(r_seed))
                j += 1
                r_seed = r_seed + j
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

        for i in range(1, len(r.reactants)+1):
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

        i = 1
        for p in products:
            f.write("variable newatoms"+str(i)+" equal floor(${counter"+str(i)+"}/2)\n")
            f.write("if '${counter"+str(i)+"} > 0' then &\n")
            f.write("'fix depositatoms all deposit ${newatoms"+str(i)+"} "+str(r_dictionary[p])+" 1 5748 region box near 2.0' \n")
            f.write("\n")
            i += 1

        #"variable newatoms equal floor(${bondcounter}/2)",
        #"if '${bondcounter} > 0' then &",
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