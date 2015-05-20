#!/bin/bash
####################################
# Prepare ligands for AutoDock Vina.
####################################
input_filename=''
output_filename=''
pH=''
split='false'
while getopts 'h:i:o:s:p:P' flag
do
  case "${flag}" in
    h) echo "usage: vina_prep.sh [-s] -i input_filename -o output_filename [-p pH] [-P]"
       echo ""
       echo "options:"
       echo "-s     split output files into one file per molecule"
       echo "-p     pH to add hydrogens"
       echo "-P     Set pH to 7.4"
       exit
       ;;
    i) input_filename=${OPTARG} ;;
    o) output_filename=${OPTARG} ;;
    s) split='true' ;;
    p) pH=${OPTARG} ;;
    P) pH='7.4' ;;
    *) echo "Unexpected option ${flag}"; exit 1 ;;
  esac
done

if [[ -z ${input_filename} ]]
then
  echo "Must provide input filename"
  exit 1
fi
if [[ -z ${output_filename} ]]
then
  echo "Must provide output filename"
  exit 1
fi
if [[ -z ${pH} ]]
then 
  echo "Must set pH"
  exit 1
fi

opts="-c -r -p ${pH} -e 1"
if [[ ${split} = 'true' ]]
then
  opts="${opts} -m"
fi

echo obabel ${input_filename} -O ${output_filename} ${opts}
obabel ${input_filename} -O ${output_filename} ${opts}
