#!/bin/sh

rotatephi=/hadron/bartek/code/rotatephi/rotatephi

if [ -z ${1} -o -z ${2} -o -z ${3} -o -z ${4} -o -z ${5} ]; then
  echo "usage: ./dorotate.sh <idir> <odir> <logdir> <L> <T>"
  echo "A logfile will be written to <logdir>, with its name derived from the basename of <odir>"
  exit 1
fi

if [ ! -d ${1} -o ! -d ${2} -o ! -d ${3} ]; then
  echo "one of the specified directories does not seem to exist!"
  exit 2
fi

logfile=${3}/$( basename ${2} ).log

echo $logfile

date > ${logfile}

for i in ${1}/*.dat; do
  # extract scalar configuration number
  cnum=$( basename $i | awk --field-separator '.' '{print $(NF-1)}' | awk --field-separator '_' '{print $NF}' )
  ${rotatephi} -i ${i} -o ${2}/scalar.$cnum -L ${4} -T ${5} | tee -a ${logfile}
done

date >> ${logfile}
