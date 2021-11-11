#! /bin/bash

usage() { echo $1; echo ;echo "Usage: ${0: -6} [-o dump.output] [-t SymulationTime] [-r RandomSeed] in.lammps " 1>&2; 
echo "Type: ${0: -6} -h for help" 1>&2; exit 1; }

# inizialise variables
ovito=false
time_value=1500
num_atoms=5
rand=10
help=false

while getopts ":o:t:r:a:h" opt; do
    case "${opt}" in
    
        o)  ovito=true
            dump=${OPTARG}
            [[ $dump == *"dump."* ]] || usage "Error: dump file needed"
            ;;

        t)  time_value=${OPTARG}
            ;;

        r)  rand=${OPTARG}
            ;;

        a)  num_atoms=${OPTARG}
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
env OMP_NUM_THREADS=16 lmp -sf omp -in $1 -var time_value $time_value -var rnseed $rand -var num_atoms $num_atoms

# to be on the safe side
sleep 1

# create a output directory if not exist
mkdir -p Output

# remove previous dump file whit same name to avoid errors
rm -r Output/$dump 2> /dev/null
rm Output/dump.* 2> /dev/null
rm Output/log.* 2> /dev/null

# move the output files to the correct directory
mv $dump Output/    2> /dev/null
mv dump.* Output/   2> /dev/null
mv log.* Output/    2> /dev/null
 
# run ovito on the dump file
if $ovito 
then
    setsid ovito Output/$dump
fi

exit 0

