#!/bin/bash
##############################################
# Run blast against sequences in the PDB
#
# usage: blast.sh query_filename [n_jobs=1]
###############################################
if [ -z $1 ]
then
  echo "Must provide query sequence filename."
  exit 1
fi

if [ -z $2 ]
then
  n_jobs=1
else
  n_jobs=$2
fi

echo "query\tsubject\tidentity\talignment_length\tmismatches\tgap_openings\tquery_start\tquery_end\tsubject_start\tsubject_end\te_value\tbit_score\tquery_length\tsubject_length\tsubject_title"
blastp -query $1 -db pdb -outfmt '6 std qlen slen stitle' -num_threads ${n_jobs}
