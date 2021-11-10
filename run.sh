#! /bin/bash

usage() { echo $1; echo ;echo "Usage: ${0: -6} [-o dump.output] [-t SymulationTime] [-r RandomSeed] in.lammps " 1>&2; 
echo "Type: ${0: -6} -h for help" 1>&2; exit 1; }

# inizialise variables
ovito=false
set_time=false
time_value=1500
rand=10
help=false

while getopts ":o:t:r:h" opt; do
    case "${opt}" in
    
        o)  ovito=true
            dump=${OPTARG}
            [[ $dump == *"dump."* ]] || usage "Error: dump file needed"
            ;;

        t)  set_time=true
            time_value=${OPTARG}
            ;;

        r)  rand=${OPTARG}
            ;;

        h)  help=true
            ;;            
            
        \?)  usage "Error: Invalid option"
            ;;
            
        : ) usage "Error: Option '$OPTARG' requires an argument" 
            ;;
    esac
done
shift $((OPTIND-1))

# display helping manual 
if $help
then
    cat $HOME/bin/instructions.txt
    exit 0
fi

# abort command, too many arguments
if [ $# -gt 1 ]
then 
    usage "Error: Too many argumets given"
fi

# abort command, too few arguments
if [ -z "$1" ] 
then
    usage   "Error: Missing input file"
fi

# run lammps input script with apropriate symulation time and random seed
env OMP_NUM_THREADS=16 lmp -sf omp -in $1 -var time_value $time_value -var rnseed $rand

# to be on the safe side
sleep 1

# run ovito on the dump file
if $ovito 
then
    setsid ovito $dump
fi

exit 0

