#!/bin/bash

OUTFILE="p_vs_alpha0.txt"

# Print Header:
echo "alpha0; p_BZ10_q0.5; p_BZ10_q0.025; p_BZ10_q0.975; p_BZ100_q0.5; p_BZ100_q0.025; p_BZ100_q0.975; p_KZ_q0.5; p_KZ_q0.025; p_KZ_q0.975" > ${OUTFILE}

# Loop on log_*.out files:
for FILE in $(ls log_*.out);
do
    PVALUE=$(grep "alpha0" ${FILE} | gawk 'BEGIN{FS=" "} {print $4}')
    gawk -v p=${PVALUE} '
        BEGIN{
            FS=" "; CNT=0
        } 
        $1=="P-value:" {
         gsub(",","",$4); 
         gsub(",","",$7); 
         gsub("\\[","",$7); 
         gsub("\\]","",$8);
	 values[(CNT * 3) + 1] = $4;
	 values[(CNT * 3) + 2] = $7;
	 values[(CNT * 3) + 3] = $8;
         CNT = CNT + 1;
        }
        END{
            print p";"values[1]"; "values[2]"; "values[3]"; "values[4]"; "values[5]"; "values[6]"; "values[7]"; "values[8]"; "values[9]
        }' ${FILE} | tr '\n' '; ' >> ${OUTFILE}
    echo "" >> ${OUTFILE}
done
