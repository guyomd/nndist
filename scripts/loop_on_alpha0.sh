#!/bin/bash

CONFIG_TEMPLATE=config_cartenat.txt
NNDIST=../nndist/bin/nndistance
NNDECLUST=../nndist/bin/nndecluster

for ALPHA in -1.0 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0;
do

    echo -e "\n\n##### RUN TEST FOR ALPHA = ${ALPHA} #####\n"

    ALSTR=$( echo ${ALPHA} | tr . p)
    sed s/"parameter_alpha0: 0.0"/"parameter_alpha0: ${ALPHA}"/ ${CONFIG_TEMPLATE} > config.txt
    sed -i s/"output_file: output.txt"/"output_file: output_${ALSTR}.txt"/ config.txt

    ${NNDIST} config.txt > log_${ALSTR}.out
    ${NNDECLUST} config.txt >> log_${ALSTR}.out
done
rm -f config.txt
