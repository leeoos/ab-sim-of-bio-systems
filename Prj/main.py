from tkinter import Tk     
from tkinter.filedialog import askopenfilename

#       -- Basic GUI --
def get_file():
    Tk().withdraw() 
    filename = askopenfilename()
    return filename


#       -- Create LAMMPS input  --
import os
import os.path
import time
from libsbml import *

def read_sbml():

    print("\nSelect a SBLM Model")
    filename = get_file()

    print("\nReading SBML file...")
    time.sleep(1)
    document = readSBML(filename)
    
    if document.getNumErrors() > 0:
        print("Encountered the following SBML errors:" )
        document.printErrors()
        return 1

    level = document.getLevel()
    version = document.getVersion()
    short_filename = filename[filename.rfind('/')+1:]

    print("\n"
                            + "File: " + short_filename
                            + " (Level " + str(level) + ", version " + str(version) + ")" )

    model = document.getModel()

    if model.isSetSBOTerm():
        print("      model sboTerm: " + model.getSBOTerm() )
        time.sleep(1)

    return model

def make_lmp():
    
<<<<<<< HEAD
        # set time steps 
        timestep 0.01   # seconds
        
        # This command sets parameters that affect
        # the building of pairwise neighbor lists
        neighbor 0.001 bin
        neigh_modify every 10 delay 100
        
        # print thermodinamic inf every N timesteps
        thermo 100
        
        # fix ID group-ID bond/create Nevery itype jtype Rmin bondtype keyword values
        # this fix will attempt to create new bond btw atoms of 
        # type 1 and 2 every Nevery timestep
        fix bonds all bond/create 10 1 2 1.0 1 prob 0.5 85784
        
        # set velocity for all atoms
        velocity all create 300.0 4928459 rot yes dist gaussian 
        
        # perform plain time integration 
        # to update position and velocity
        # and simulate Brownioan motion
        fix 1 all nve\n
        fix 2 all langevin 300.0 300.0 10.0 904297  
        
        # compute if atoms has a bonds
        # and total number of bonds btw all atoms
        compute 1 agents property/atom nbonds
        compute 2 agents reduce sum c_1            
        thermo_style custom step temp pe c_2
        run 0
        # this lines are necessary to insure that the “hasbond” and "newatoms" 
        # variables are current when the group command invokes it.
        
        # hasbond : boolean = true if atom I has a bond with atom J
        variable hasbond atom "c_1 > 0.0"
        
        # bondcounter : int = N total number of bonds in the sim
        variable bondcounter equal ceil(c_2) 
        
        # print themo info every timestep 
        thermo_style custom step temp pe v_bondcounter
        
        # dumps atoms information 
        dump 1 all custom 10 dump.out id x y z type 
        
        """
    )
    f_in.write(
        """
        #       --- LOOP---
        
        label loop\n
        variable step loop ${duration}   # loop length
        
        # create new atoms only if new bonds have been made 
        # the num of new atoms is linked to the number of new bonds as follow:
        variable newatoms equal floor(${bondcounter}/2)
        if "${bondcounter} > 0" then &
        "fix depositatoms all deposit ${newatoms} 3 1 5748 region box near 2.0" 
        # fix ID group-ID deposit N type M seed keyword values
        
        # assing all atoms that have a bond to the garbage group
        group garbage dynamic all every 1 var hasbond
        
        # append new values on dump file
        dump_modify 1 append yes
        
        # perform n steps in loop
        run 100
=======
    os.system('mkdir -p LAMMPS')

    # rand seed starter
    r_seed = 5783

    model = read_sbml()

    if model == 1:
        os._exit(1)

    print("\nCreating lammps file as in.lmp in LAMMPS/ dir...")
    time.sleep(1)

    # info
    species = model.getNumSpecies()

    with open('LAMMPS/in.lmp', 'w') as f:
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
        "create_box  3 box  bond/types 1 extra/bond/per/atom 10\n",

        "# create simulation walls",
        "fix xwalls all wall/reflect xlo EDGE xhi EDGE",
        "fix ywalls all wall/reflect ylo EDGE yhi EDGE",
        "fix zwalls all wall/reflect zlo EDGE zhi EDGE\n"]

        f.writelines(["%s\n" % item  for item in sim_box])

        # AGENTS PROPRETIES AND FORCE FIELDS
        f.write("\n#       --- AGENTS PROPRETIES AND FORCE FIELDS ---\n")
        f.write("# creation of atoms of types in randoms spots inside the box\n")
        
        for i in range(1, species+1):
            f.write("create_atoms" + "    " + str(i) + " random ${atoms} " + str(r_seed+i) + " box")
            f.write("\n")

        ag_prop =[
        "\n# atoms mass",
        "mass 1 10.948",
        "mass 2 10.467",
        "mass 3 10.578   # new atoms generated by type 1-2 bond\n",
        
        "# assing atoms to cerian groups",
        "group g1 type 1 ",
        "group g2 type 2",
        "group agents  union g1 g2\n",
        
        "# force fields style and coefficient",
        "pair_style zero 5.0",
        "pair_coeff * *\n",
        
        "# bond style and coefficients",
        "bond_style  harmonic",
        "bond_coeff * 100 1.1\n"]

        f.writelines(["%s\n" % item  for item in ag_prop])
        
        # SIMULATION 
        sim = [
        " \n#      --- SIMULATION ---\n",

        "# set time steps ",
        "timestep 0.01   # seconds\n",
        
        "# This command sets parameters that affect",
        "# the building of pairwise neighbor lists",
        "neighbor 0.001 bin",
        "neigh_modify every 10 delay 100\n",
        
        "# print thermodinamic inf every N timesteps",
        "thermo 100 \n",
        
        "# fix ID group-ID bond/create Nevery itype jtype Rmin bondtype keyword values",
        "# this fix will attempt to create new bond btw atoms of ",
        "# type 1 and 2 every Nevery timestep",
        "fix bonds all bond/create 10 1 2 1.0 1 prob 0.5 85784 \n",
        
        "# set velocity for all atoms",
        "velocity all create 300.0 4928459 rot yes dist gaussian \n",
        
        "# perform plain time integration ",
        "# to update position and velocity",
        "# and simulate Brownioan motion",
        "fix 1 all nve",
        "fix 2 all langevin 300.0 300.0 10.0 904297\n",
        
        "# compute if atoms has a bonds",
        "# and total number of bonds btw all atoms",
        "compute 1 agents property/atom nbonds",
        "compute 2 agents reduce sum c_1 ",           
        "thermo_style custom step temp pe c_2",
        "run 0",
        "# this lines are necessary to insure that the “hasbond” and 'newatoms' ",
        "# variables are current when the group command invokes it.\n",
        
        "# hasbond : boolean = true if atom I has a bond with atom J",
        "variable hasbond atom 'c_1 > 0.0'\n",
        
        "# bondcounter : int = N total number of bonds in the sim",
        "variable bondcounter equal ceil(c_2) \n",
        
        "# print themo info every timestep ",
        "thermo_style custom step temp pe v_bondcounter\n",
        
        "# dumps atoms information ",
        "dump 1 all custom 10 dump.out id x y z type \n"]

        f.writelines(["%s\n" % item  for item in sim])

        # LOOP 
        loop = [
        "\n#       --- LOOP---\n",
        
        "label loop",
        "variable step loop ${duration}   # loop length\n",
        
        "# create new atoms only if new bonds have been made", 
        "# the num of new atoms is linked to the number of new bonds as follow:",
        "variable newatoms equal floor(${bondcounter}/2)",
        "if '${bondcounter} > 0' then &",
        "'fix depositatoms all deposit ${newatoms} 3 1 5748 region box near 2.0' ",
        "# fix ID group-ID deposit N type M seed keyword values\n",
        
        "# assing all atoms that have a bond to the garbage group",
        "group garbage dynamic all every 1 var hasbond\n",
        
        "# append new values on dump file",
        "dump_modify 1 append yes\n",
        
        "# perform n steps in loop",
        "run 100\n",
>>>>>>> python
    
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

        f.writelines(["%s\n" % item  for item in loop])

        f.flush()
        os.fsync(f)

import keyboard

def run_():

    str_in = 'start'
    while True:
        os.system('clear')
        if (str_in == 'start'):
            print("\n   Welcome: \n")
            print("Press 1 to generate a input script\n")
            print("Type q to exit\n")
            str_in = input("> ")
            if (str_in == '1') : 
                make_lmp()
                str_in = '0'

        if (str_in != '' and str_in != 'q'):
            os.system('clear')
            print("\n   Action Menu: \n")
            print("Press 1 to look at the lammps script\n")
            print("Press 2 to run the script\n")
            print("Press 3 to restart\n")
            
            print("Type q to exit\n")

        if (str_in != 'q'): str_in = input("> ")

        
        if (str_in == '1') : 
            os.system('less ./LAMMPS/in.lmp')
            str_in = '0'

        elif (str_in == '2') :
            os.system('./run.sh -o dump.out LAMMPS/in.lmp')
            str_in = '0'

        elif (str_in == '3'):
            os.system('python3 main.py')
            return 0
    
        elif (str_in == 'q'):
            print("\nexiting...\n")
            time.sleep(1)
            os.system('clear')
            return 0
        
        elif (str_in == '') : None

        else: print("Error\n")


if __name__ == '__main__':
    run_()
