#!/bin/bash
####################################
# Prepare ligands for AutoDock Vina.
####################################
input_filename=''
split='false'
while getopts 'hi:s' flag
do
  case "${flag}" in
    h) echo "usage: vina_prep.sh [-s] -i input_filename"
       echo ""
       echo "options:"
       echo "-s     split output files into one file per molecule"
       exit
       ;;
    i) input_filename=${OPTARG} ;;
    r) receptor='true' ;;
    s) split='true' ;;
    *) echo "Unexpected option ${flag}"; exit 1 ;;
  esac
done

if [[ -z ${input_filename} ]]
then
  echo "Must provide input filename"
  exit 1
fi

opts=''
if [[ ${split} = 'true' ]]
then
  opts="${opts} -m"
fi

output_basename=$(basename ${input_filename} | sed 's/\..*$//g')

echo obabel ${input_filename} -O ${output_basename}.pdbqt ${opts}
obabel ${input_filename} -O ${output_basename}.pdbqt ${opts}
