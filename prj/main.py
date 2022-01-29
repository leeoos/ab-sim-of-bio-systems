# /bin/python3

from fileinput import filename
from io import SEEK_CUR
from tkinter import Radiobutton, Tk     
from tkinter.filedialog import askopenfilename
from webbrowser import get

#       -- Basic GUI --
import os
import os.path
import time

def get_file():
    Tk().withdraw() 
    print("\nSelect a SBLM Model")
    filename = askopenfilename()
    if filename == () : os._exit(3)
    return filename

from lib.create import *

def run_():

    lmp_dir = 'simulation'

    filename = get_file()
    short_filename = filename[filename.rfind('/')+1:-4]

    PATH = os.path.dirname(os.path.realpath(__file__))
    os.system('mkdir -p '+PATH+'/'+lmp_dir+'/')

    str_in = 'start'

    while True:
        os.system('clear')

        if (str_in == 'start'):
            print("\n   Welcome: \n")
            print("Press 1 to generate a input script\n")
            print("Type q to exit\n")
            str_in = input("> ")

            if (str_in == '1'):
                print("\nThe actual rand seed is: 5783 \n")
                print("Type a new number if you want to change it\n")
                print("Press Enter to continue with the default\n")
                n_seed = input("> ")

                print("\nThe actual "+lmp_dir+" time is: 5000 \n")
                print("Type a new number if you want to change it\n")
                print("Press Enter to continue with the default\n")
                s_time = input("> ")
                if (s_time == "" ) : s_time = "5000"

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
            print("Press 3 to restart\n")
            
            print("Type q to exit\n")

        if (str_in != 'q'): str_in = input("> ")

        
        if (str_in == '1') : 
            os.system('less '+PATH+'/'+lmp_dir+'/in.lmp')
            str_in = '0'

        elif (str_in == '2') :
            os.system(PATH+'/run.sh -t '+s_time+' -o dump.'+short_filename+' '+PATH+'/'+lmp_dir+'/in.lmp '+ short_filename)
            str_in = '0'

        elif (str_in == '3'):
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
 

