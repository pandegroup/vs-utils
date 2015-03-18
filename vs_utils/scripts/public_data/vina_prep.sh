#!/bin/bash
##################################################
# Prepare ligands and receptors for AutoDock Vina.
##################################################
input_filename=''
output_filename=''
receptor='false'
while getopts 'hi:o:r' flag
do
  case "${flag}" in
    h) echo "usage: vina_prep.sh [-r] -i input_filename -o output_filename"
       echo ""
       echo "options:"
       echo "-r     treat molecules as receptors (rigid)"
       exit
       ;;
    i) input_filename=${OPTARG} ;;
    o) output_filename=${OPTARG} ;;
    r) receptor='true' ;;
    *) echo "Unexpected option ${flag}"; exit 1 ;;
  esac
done
if [[ ${receptor} = 'true' ]]
then
  opts='-xr'
else
  opts=''
fi
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
obabel ${input_filename} -O ${output_filename} ${opts}
