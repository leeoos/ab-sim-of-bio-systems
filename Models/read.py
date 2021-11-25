import sys
import time
import os
import os.path
from libsbml import *

def main (args):
    """Usage: readSBML filename
    """

    if len(args) != 2:
        print("Usage: readSBML filename")
        return 1

    filename = args[1]
    current = time.process_time()
    document = readSBML(filename)

    errors = document.getNumErrors()

    print()
    print("            filename: " + filename)
    print("           file size: " + str(os.stat(filename).st_size))
    print("      read time (ms): " + str(time.process_time() - current))
    print(" validation error(s): " + str(errors))
    print()
    document.printErrors()

    return errors


if __name__ == '__main__':
    main(sys.argv) 
     