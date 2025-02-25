#!/bin/bash

ABACUS_PATH=$(awk -F "=" '$1=="ABACUS_PATH"{print $2}' ../../SETENV)
ABACUS_NPROCS=$(awk -F "=" '$1=="ABACUS_NPROCS"{print $2}' ../../SETENV)
ABACUS_THREADS=$(awk -F "=" '$1=="ABACUS_THREADS"{print $2}' ../../SETENV)

for i in $( seq 0 7 )
do
    cp INPUT_$i INPUT
    OMP_NUM_THREADS=${ABACUS_THREADS} mpirun -np ${ABACUS_NPROCS} ${ABACUS_PATH} | tee md$i.output
    cp OUT.ABACUS/running_md.log running_md$i.log
done

rm INPUT

for i in $( seq 0 7 )
do
    if [[ ! -f md$i.output ]] || 
       [[ ! -f OUT.ABACUS/running_md$i.log ]] ||
       [[ ! ( "$(tail -1 OUT.ABACUS/running_md$i.log)" == " Total  Time  :"* ) ]]
    then
    	echo "job is failed!"
    	exit 1
    fi
done

echo "job is successed!"
exit 0