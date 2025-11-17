#!/bin/bash

CONFIG_TEMPLATE=config_ZBZ2020_alpha0.1.txt
NNDIST=bin/nndistance
NNDECLUST=bin/nndecluster

for ALPHA in -1.0 -0.8 -0.6 -0.4 -0.2 0.0 0.2 0.4 0.6 0.8 1.0 0.1 0.3 0.5;
do

    echo -e "\n\n##### RUN TEST FOR ALPHA = ${ALPHA} #####\n"

    sed s/"parameter_alpha0: 0.1"/"parameter_alpha0: ${ALPHA}"/ ${CONFIG_TEMPLATE} > config.txt
    sed -i s/"output_file: output_ZBZ2020_alpha0.1.txt"/"output_file: output_test.txt"/ config.txt

    ${NNDIST} config.txt > log_al_${ALPHA}.out
    ${NNDECLUST} config.txt >> log_al_${ALPHA}.out
done
rm -f config.txt
