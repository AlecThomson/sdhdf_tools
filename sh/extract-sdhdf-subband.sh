#!/bin/bash
#
# extract-sdhdf-subband.sh
#
# Program to extract a sub-band and associated data and meta-data
# from an SDHDF file and write out to a new SDHDF file
#

SCRIPT=`basename ${0}`
SB=${1}
INFILE=${2}
OUTFILE=${3}

if [ ${#} -ne 3 ]; then
  echo "Program: ${SCRIPT}"
  echo "Usage: ./${SCRIPT} <sub-band to extract [0-25]>"
  echo "                          <input_filename>"
  echo "                          <output_filename>"
  exit 0
fi

# Sanity checks -----------------------------
if [ ! -e ${INFILE} ]; then
  echo "${SCRIPT}:: ERROR: File '${INFILE}' does not exist. Exiting."
  exit 1
fi

if [ -e ${OUTFILE} ]; then
  echo "${SCRIPT}:: ERROR: File '${OUTFILE}' already exists. Exiting."
  exit 1
fi

if [ ${SB} -le 9 ]; then
  SB="0${SB}"
elif [ ${SB} -gt 25 ]; then
  echo "${SCRIPT}:: ERROR: Max sub-band is 25"
  exit 1
fi

h5_objs=`h5ls ${INFILE} | awk '{print$1}' | grep -v SB`
sb_objs="band_SB${SB}"

# Main -------------------------------------
echo -e "\n${SCRIPT}:: Processing sub-band '${SB}'...\n"

for obj in ${h5_objs} ${sb_objs}
do
  echo "${SCRIPT}:: Extracting ${obj}..."
  h5copy -i ${INFILE} -o ${OUTFILE} -s "/${obj}" -d "/${obj}"
  stat="${?}"
  if [ ${stat} -ne 0 ]; then
    echo "${SCRIPT}:: Extracting from file '${INFILE}' failed. Exiting."
    exit 1
  fi
done 

if [ -e ${OUTFILE} ]; then
  echo -e "\n${SCRIPT}:: New file '${OUTFILE}' contains the following:\n"
  h5ls ${OUTFILE}
else
  echo "${SCRIPT}:: Extracting from file '${INFILE}' failed. Exiting."
fi

