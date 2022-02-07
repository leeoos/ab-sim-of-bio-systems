# /bin/python3

from io import SEEK_CUR
from tkinter import Radiobutton, Tk     
from tkinter.filedialog import askopenfilename

#       -- Basic GUI --
import os
import os.path
import time

def get_file():
    Tk().withdraw() 
    os.system('clear')
    print("\nSelect a SBLM Model")
    filename = askopenfilename()
    if filename == () : os._exit(3)
    return filename

# importing create script to traduce sbml into lammps
from lib.create import *

# -- to implement as a GUI
def run_():

    filename = get_file()
    short_filename = filename[filename.rfind('/')+1:-4]
    lmp_dir = 'simulation/'+short_filename

    PATH = os.path.dirname(os.path.realpath(__file__))
    os.system('mkdir -p '+PATH+'/'+lmp_dir+'/')

    str_in = 'start'

    while True:
        os.system('clear')

        if (str_in == 'start'):
            print("\n   Welcome: \n")
            print("Press 1 to generate a input script\n")
            print("Type q to exit\n")
            str_in = input("> ").lower()
            os.system('clear')

            if (str_in == '1'):
                print("\nThe actual rand seed is: 5783 \n")
                print("Type a new number if you want to change it\n")
                print("Press Enter to continue with the default\n")
                n_seed = input("> ")
                os.system('clear')

                print("\nThe current simulation time is set to 500 steps \n")
                print("Type a new number if you want to change it\n")
                print("Press Enter to continue with the default\n")
                s_time = input("> ")
                os.system('clear')

                if (s_time == "" ) : s_time = "500"

                if (n_seed != "") : make_lmp(r_seed=int(n_seed), 
                    lmp_file_path=PATH+'/'+lmp_dir+'/in.lmp', sbml_filename = filename)
                else : make_lmp(lmp_file_path=PATH+'/'+lmp_dir+'/in.lmp', 
                    sbml_filename = filename)
                    
                str_in = '0'


        if (str_in != '' and str_in != 'q'):
            os.system('clear')
            print("\n   Action Menu: \n")
            print("Press 1 to look at the lammps script\n")
            print("Press 2 to run the script\n")
            print("Press 3 to run the script with ovito\n")
            print("Press r to restart\n")
            print("Type q to exit\n")

        if (str_in != 'q'): str_in = input("> ").lower()
        
        if (str_in == '1') : 
            os.system('less '+PATH+'/'+lmp_dir+'/in.lmp')
            str_in = '0'

        elif (str_in == '2') :
            os.system(PATH+'/run.sh -t '+s_time+' '+PATH+'/'+lmp_dir+'/in.lmp '+ short_filename)
            str_in = '0'

        elif (str_in == '3') :
            os.system(PATH+'/run.sh -t '+s_time+' -o dump.'+short_filename+'.out '+PATH+'/'+lmp_dir+'/in.lmp '+ short_filename)
            str_in = '0'

        elif (str_in == 'r'):
            os.system('python3 '+ PATH +'/main.py')
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
 

