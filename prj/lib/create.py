# /bin/python3
# This scrip create a LAMMPS input file 
# based on the SBM model given in input

#       --- CLASSES ---
from functools import reduce
from decimal import Decimal
from mimetypes import init
from random import randint
import math
import os

# get info about the
class CompartmentClass:

    def __init__(self,model) -> None:
        self.num_comp = model.getNumCompartments()
        if (self.num_comp) > 1 :
            print(  "\nWarning: more than one compartment detected. "
                    +"This may affect the correct execution of the simulation. "
                    +"Consult AB-Sim-Of-Bio-Systems/resources/Tesi.pdf for more info "
                )
        self.csize_dict = self.__make_ccsize_dict(model)
        pass

    def __make_ccsize_dict(self, model):
        csize_dict = {}
        for i in range(self.num_comp):
            comp = model.getCompartment(i)
            csize = comp.getSize()
            if (math.isnan(csize)) : csize = 1
            csize_dict[comp.getId()] = csize
        return csize_dict

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
            csize = C.csize_dict[specie.getCompartment()] 
            if (math.isnan(value)) : scale = False; break
            else : initial_atoms.append(math.ceil(value* pow(10,-(self.__fexp(value)))) * int(csize))
        if not(scale):
            initial_atoms = []
            for i in range(self.num_of_species):
                specie = model.getSpecies(i)
                value = specie.getInitialAmount()
                if (math.isnan(value)) : 
                    print('Error: no initial ammount/concentration given')
                    os._exit(3)
                else :
                    initial_atoms.append(value)     ## CHECK!!!
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
            growth_rate = randint(100, 500)
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
        test = 'Tests/simple.xml'
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
    
        ###     SET UP OF INPUT VARIABLES       ###
        set_up =[
        "\n#       --- SET UP OF INPUT VARIABLES ---\n",

        "# steps : int = number of steps to execute",
        "variable time_value index 10000  # default value: 10000",
        "variable steps equal ${time_value}\n",

        "# duration : int = simulation steps * dump interval",
        "#for the current run",
        "variable duration equal ${steps}*0.1\n",

        "# uncomment the last 2 lines of this section and",
        "# substitue the number of atoms in the create_atoms command with ${atoms} ",
        "# to modify the initial amount of atoms in the simulation;",
        "# prj/run.sh -h for more info",
        "# Note: doing so the initial amount will be the same for all atoms types",
        "# atoms : int = N number of atoms of each type to generate ",
        "#variable num_atoms index 5",
        "#variable atoms equal ${num_atoms}\n"]
        f.writelines(["%s\n" % item  for item in set_up])

        if (S.total_initial_atoms <= 1 ) : expected_types = S.num_of_species + 1
        else : expected_types = S.num_of_species
    
        ###      SIMULATION BOX PROPERTIES      ###
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
        "create_box  "+str(expected_types)+" box  bond/types 1 extra/bond/per/atom 100\n", # ATTENZIONE ALL'EXTRA BOND TYPE

        "# create simulation walls",
        "fix xwalls all wall/reflect xlo EDGE xhi EDGE",
        "fix ywalls all wall/reflect ylo EDGE yhi EDGE",
        "fix zwalls all wall/reflect zlo EDGE zhi EDGE\n"]
        f.writelines(["%s\n" % item  for item in sim_box])

        ###      AGENTS PROPRETIES AND FORCE FIELDS     ###
        f.write("\n#       --- AGENTS PROPRETIES AND FORCE FIELDS ---\n\n")
        f.write("# creation of atoms of types in randoms spots inside the box\n")
        
        for k in list(S.dictionary.keys()):
            types = ""
            if k in R.reactants : types = types + "     | reactant |"
            if k in R.products : types = types + "      | product  |"
            f.write("create_atoms"+"    "+ str(S.dictionary[k][0]) +" random "+ str(S.dictionary[k][2]) +
                        " "+ str(r_seed) + " box  # "+ k + types)
                        
            r_seed = r_seed + 1
            f.write("\n")
        
        # in case there are 1 or less atom at the beginning
        # create ghost atoms to inizialize velocity
        if (S.total_initial_atoms <= 1): 
            f.write("create_atoms"+"    "+ str(S.num_of_species+1) +" random 2 "
                        + str(r_seed) + " box  # ghost atom to inizialize velocity")
            r_seed = r_seed + 1
            f.write("\n")
        else : pass

        f.write("\n# atoms mass\n")
        for i in list(S.dictionary.values()) :
            f.write("mass " + str(i[0]) + " 1.0")
            f.write("\n")
        
        # assign mass to ghost atoms
        if (S.total_initial_atoms <= 1 ) : 
            f.write("mass " + str(S.num_of_species+1) + " 0.5   # ghost atom")
            f.write("\n")
        else : pass

        # make a set of perishible atoms
        perishable = set()
        for r_id in list(R.reactions.keys()): 
            if (len(R.reactions[r_id][0]) ==  1) :
                key = R.reactions[r_id][0][0]
                perishable.add(S.dictionary[key][0])
            else : pass

        # create a grup to dump in which ghost atoms are not incuded
        # and one (ghost) in which are, so they can be easly delated
        if (S.total_initial_atoms <= 1) :
            f.write(
                "\n# create a grup to dump in which ghost atoms are not incuded \n" +
                "# and one in which are, so they can be easly delated\n"
            )
            f.write("group to_dump empty\n")
            f.write("group ghost type "+ str(S.num_of_species+1) +"\n")
        else : pass
 
        ag_prop = [
        "\n# force fields style and coefficient",
        "pair_style zero 5.0",
        "pair_coeff * *\n",
        
        "# bond style and coefficients",
        "bond_style  harmonic",
        "bond_coeff * 100 1.1\n"]
        f.writelines(["%s\n" % item  for item in ag_prop])
        
        ###      SIMULATION         ###
        sim1 = [
        " \n#      --- SIMULATION ---\n",

        "# set the time step for the simulation",
        "# it is directly proportional to the atoms velocity",
        "# AB-Sim-Of-Bio-Systems/resources/Tesi.pdf for more info",
        "timestep 0.001   # abstract units : seconds\n",
        
        "# this command sets parameters that affect",
        "# the building of pairwise neighbor lists",
        "neighbor 0.001 bin",
        "neigh_modify every 10 delay 100\n",
        
        "# print thermodinamic inf every N timesteps",
        "thermo 100 \n"]
        f.writelines(["%s\n" % item  for item in sim1])
        
        if (len(R.combinations) > 0) : 
            f.write("# this fix will attempt to create new bond btw atoms of type i and j every N timestep\n")
            f.write("# fix ID group-ID bond/create Nevery itype jtype Rmin bondtype keyword values\n")
        else : pass
        
        # TODO : EDIT TYPE OF BOND AND PROB 
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
        "fix kin all nve",
        "fix lgv all langevin 300.0 300.0 10.0 904297\n",

        ###     COMPUTES        ###
        "\n#        ---COMPUTES--- "]
        f.writelines(["%s\n" % item  for item in sim2])

        # compute perishable atoms
        if not(perishable == set()) :
            to_write = (
                "\n# comupe which atoms should be inhibid \n" +
                "# by checking their type and kin energy \n" +
                "compute atype all property/atom type \n" +
                "compute akin all ke/atom \n\n" +
                "# toInhibid : Boolean = true if atom I is\n" +
                "# perishable and has a low kin value \n"
            )
            f.write(to_write)
            toInhibid = "variable toInhibid atom '("
            for i in perishable : 
                toInhibid = toInhibid + " c_atype == "+ str(i) + " || "
            toInhibid = toInhibid[:-4] + ") & \n&& c_akin < 100.0'" 
            f.write(toInhibid + "\n")
        else : pass

        # compute ghost atoms
        to_dump = "all"
        if (S.total_initial_atoms <= 1) :
            if (perishable == set()) :
                to_write = (
                    "\n# comupte which atoms are ghost atoms \n"+
                    "# by checking their type \n"+
                    "compute atype all property/atom type \n"
                )
                f.write(to_write)
            to_dump = "to_dump"
            to_write = (
                "\n# isToDump : boolean = true if atom I type is between 1 and " +
                str(S.num_of_species) + 
                "\nvariable isToDump atom 'c_atype < "+ str(S.num_of_species+1) +"'\n" 
            )
            f.write(to_write)
        else: 
            to_dump = "all"
        
        # compute interactions
        if (len(R.combinations) > 0) :
            if ((S.total_initial_atoms >= 1) and (perishable == set())) :
                f.write("\n# compute the atom type for each atom \n")
                f.write("compute atype all property/atom type \n")

            f.write("\n# compute number of bonds for each atoms \n")
            f.write("compute nbond all property/atom nbonds \n\n")

            # check for atoms that has to be delated
            f.write("# toDelate : boolean = true if the atom I has reached\n# the max num of bonds for its type")
            toDelate = "\nvariable toDelate atom '("
            next_line = 1
            for r_id in list(R.reactions.keys()):
                num_bonds = int(len(R.reactions[r_id][0])-1)
                if len(R.reactions[r_id][0]) <= 1 : pass
                else :
                    for atype in (R.reactions[r_id][0]) :
                        if (next_line % 3 != 0):
                            toDelate = (
                                toDelate + "c_atype == "+ str(S.dictionary[atype][0]) 
                                +" && " + "c_nbond >= " + str(num_bonds) + ") || ("
                            )
                        else :
                            toDelate = (
                                toDelate + "c_atype == "+ str(S.dictionary[atype][0]) 
                                +" && " + "c_nbond >= " + str(num_bonds) + ") & \n|| ("
                            )
                        next_line += 1
            toDelate = toDelate[:-4] + "' \n"
            f.write(toDelate)
            
            # count new atoms by computing the number of different bond in 
            # the simulation at the current timestep
            f.write("\n# count type of bond present in the \n")
            f.write("# simulation at the current timestep \n")
            compute_reduce = ""
            for i, r_id in enumerate(list(R.reactions.keys())):
                if len(R.reactions[r_id][0]) <= 1 : pass
                else :
                    atom_type = S.dictionary[R.reactions[r_id][0][0]][0]
                    num_bonds = len(R.reactions[r_id][0]) - 1
                    to_write = (
                        "variable "+ str(r_id) + "Bonds atom '(c_atype == "+ str(atom_type) 
                        +" && c_nbond == "+ str(num_bonds) +")' \n" +
                        "compute count"+ str(r_id) +" all reduce sum v_"+ str(r_id) + "Bonds \n\n"
                    )
                    f.write(to_write)
                    if ((i+1) % 4 == 0) : 
                        compute_reduce = compute_reduce + "c_count"+ str(r_id) +" &\n"
                    else : 
                        compute_reduce = compute_reduce + "c_count"+ str(r_id) +" "
            to_write = (
                "# this lines are necessary to ensure that the compute \n"+
                "# reduce reference is current when invoked.\n" +
                "thermo_style custom step temp pe "+ compute_reduce + "\n"+
                "run 0 \n\n"
            )
            f.write(to_write)
            
            # compute number of new produts for each reaction
            new_products = []
            for r_id in list(R.reactions.keys()):
                if len(R.reactions[r_id][0]) <= 1 : pass
                else :
                    f.write("# produts for reaction "+ str(r_id) +"\n")
                    for p in R.reactions[r_id][1] :
                        new_atom = str(S.dictionary[p][0]) +"_"+ str(S.dictionary[p][4])
                        S.dictionary[p][4] += 1
                        f.write("variable newAtoms"+ new_atom +" equal c_count" + str(r_id) + "\n")
                        new_products.append("newAtoms"+ new_atom)
                    f.write("\n")
        else : 
            f.write("\n")
            pass
    
        ###      DINAMICS GROUPS       ###
        f.write("\n#        ---DINAMICS GROUPS--- \n")
        f.write("# note: this section could be empty\n")
       
        # update garbage dynamic group
        if (len(R.combinations) > 0) :
            to_write = (
                "\n# assing all atoms that have shold be delated to  \n"+
                "# the garbage group and update it every N timesteps \n"+
                "group garbage dynamic all every 10 var toDelate \n" 
            )
            delate_garbage = "delete_atoms group garbage bond yes mol yes compress no"
            f.write(to_write)
        else : delate_garbage = ""
       
        if not(perishable == set()) :
            to_write = (
                "\n# assing some perishable atoms to the inhibit group \n"+
                "# update group every N timesteps \n"+ 
                "group inhibid dynamic all every 10 var toInhibid \n"
            )
            delate_perish = "delete_atoms group inhibid compress no"
            f.write(to_write)
        else : delate_perish = ""

        # dump all atoms exept ghost atoms
        if(S.total_initial_atoms <= 1):
            to_write = (
                "\n# assing all atoms of the correct type to the 'to-dump' group\n"+
                "group to_dump dynamic all every 10 var isToDump \n"
            )
            f.write(to_write)
        else: pass

        ###      DUMP FILE       ###
        f.write("\n\n#        ---DUMP FILE--- \n\n")

        # make a dump file to sote sim info every 10 timestep 
        f.write("# dump atoms information every 10 timestep\n")
        f.write("# dump ID group-ID style N file args\n")
        f.write("dump d1 "+ to_dump + " custom 10 dump."+ short_filename +".out id x y z type \n")

        ###      DEPOSIT       ### 
        f.write("\n\n#        ---DEPOSIT--- \n\n")
        to_write = (
            "# insert a single atom or molecule into the simulation domain \n" +
            "# every N timesteps until M atoms or molecules have been inserted \n"+
            "# fix ID group-ID deposit M type N seed keyword values \n"
        )
        f.write(to_write)
        for r_id in list(R.reactions.keys()):
            if (R.reactions[r_id][1] == []) :  pass # reactions with no produts
            else :
                for p in R.reactions[r_id][1] :
                    atom = str(S.dictionary[p][0]) +"_"+ str(S.dictionary[p][4])
                    S.dictionary[p][4] += 1
                    p_id = str(S.dictionary[p][0])
                    growth = S.dictionary[p][3]
                    if (len(R.reactions[r_id][0]) <= 1) :
                        f.write("fix deposit"+ atom +" all deposit ${duration} "+ p_id +" "
                                    + str(growth) +" "+str(r_seed) +" region box near 2.0\n")
                        r_seed += 1
                    else : pass

        ###      RUN       ### 
        f.write("\n\n#        ---RUN--- \n")

        if (len(R.combinations) > 0) :
            to_write = (
                "\n# uncomment this line to append\n"+
                "# new values on dump file in case \n"+
                "# it is not updating properly \n"+
                "# dump_modify d1 append yes \n\n"+
                "# allow atoms lost during the simulation \n"+
                "thermo_modify lost ignore \n\n"+
                "run ${steps} every 100 &\n"+
                "   'variable position equal ceil(random(1,1000000000,10)/100000.0) ' & \n" 
            )
            f.write(to_write)

            # TODO : number of atom to create
            create = [] ; pos = 0
            for r_id in list(R.reactions.keys()):
                if (R.reactions[r_id][1] == []) :  pass # reactions with no produts
                else :
                    for p in R.reactions[r_id][1] :
                        if (len(R.reactions[r_id][0]) <= 1) : pass
                        else:
                            create.append("   'create_atoms "+ str(S.dictionary[p][0]) 
                                            +" random ${"+ new_products[pos] +"} ${position}+"
                                            + str(pos) +" box ' &"
                                        )
                            pos += 1
            f.writelines(["%s\n" % item  for item in create])

            to_write = (
                "   '"+ delate_garbage +"' & \n" +
                "   '"+ delate_perish +"' \n"
            )
            f.write(to_write)
        elif not(perishable == set() ) :
            to_write = (
                "# uncomment this line to append\n"+
                "# new values on dump file in case \n"+
                "# it is not updating properly \n"+
                "# dump_modify d1 append yes \n\n"+
                "# allow atoms lost during the simulation \n"+
                "thermo_modify lost ignore \n"+
                "run ${steps} every 100 &\n"+
                "   '"+ delate_perish +"' \n" 
            )
            f.write(to_write)
        else :
            f.write("run ${steps}")

        ###     END     ###
        f.write("\n#        ---END--- \n\n") 
        to_write = (
            "# print some useful informations \n"+
            "print '' \n"+
            "print 'INFO' \n" +
            "print 'Starting Atoms: "+ str(S.total_initial_atoms) +" ' \n" +
            "print 'Duration: ${duration}' \n" +
            "print 'ALL DONE' \n"
        )
        f.write(to_write)    
        
        f.flush()
        os.fsync(f)

        print("\nALL DONE\n")

import logging
import sys
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

if __name__ == '__main__':
    make_lmp()