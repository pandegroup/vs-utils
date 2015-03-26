#!/bin/bash
##################################################
# Prepare ligands and receptors for AutoDock Vina.
##################################################
input_filename=''
receptor='false'
split='false'
while getopts 'hi:rs' flag
do
  case "${flag}" in
    h) echo "usage: vina_prep.sh [-rs] -i input_filename"
       echo ""
       echo "options:"
       echo "-r     treat input as a receptor (rigid)"
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
if [[ ${receptor} = 'true' ]]
then
  opts='-xr'
fi

if [[ ${split} = 'true' ]]
then
  opts="${opts} -m"
fi

output_basename=$(basename ${input_filename} | sed 's/\..*$//g')

echo obabel ${input_filename} -O ${output_basename}.pdbqt ${opts}
obabel ${input_filename} -O ${output_basename}.pdbqt ${opts}