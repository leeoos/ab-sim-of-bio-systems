from io import SEEK_CUR
from tkinter import Radiobutton, Tk     
from tkinter.filedialog import askopenfilename

#       -- Basic GUI --
def get_file():
    Tk().withdraw() 
    print("\nSelect a SBLM Model")
    filename = askopenfilename()
    return filename

import os
import os.path
import time

from lib.create import *

def run_():

    os.system('mkdir -p models/lammps/')

    # default rand seed STARTER -> it will be incremented each time
    r_seed = 5783

    str_in = 'start'

    while True:
        os.system('clear')

        if (str_in == 'start'):
            print("\n   Welcome: \n")
            print("Press 1 to generate a input script\n")
            print("Type q to exit\n")
            str_in = input("> ")

            if (str_in == '1') :
                print("\nThe actual rand seed is: " + str(r_seed) + '\n')
                print("Type a new number if you want to change the rand seed\n")
                print("Press Enter to continue with the default\n")
                n_seed = input("> ")

                if (n_seed != '') : make_lmp(r_seed=int(n_seed), 
                    lmp_file_path='./models/lammps/in.lmp', sbml_filename = get_file())
                else : make_lmp(lmp_file_path='./models/lammps/in.lmp', 
                    sbml_filename = get_file())
                    
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
            os.system('less ./models/lammps/in.lmp')
            str_in = '0'

        elif (str_in == '2') :

            os.system('./run.sh -t 50000 -o dump.out ./models/lammps/in.lmp')

            os.system('./run.sh -o dump.out ./models/lammps/in.lmp')

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
