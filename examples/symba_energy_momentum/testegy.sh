#!/bin/bash
cases=("merger" "purehitandrun" "hitandrun" "disruption" "supercatastrophic")
for i in "${cases[@]}"
do
   ./swiftest_symba "param.${i}.in" | tail -n5 | head -n1 | awk '{print $4,$5,$6,$7,$8,$9}' | sed 's|^|'${i}'\t|'
done
