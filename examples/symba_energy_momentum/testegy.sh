#!/bin/bash
cases=("merger" "purehitandrun" "hitandrun" "disruption" "supercatastrophic")
for i in "${cases[@]}"
do
   ./swiftest_symba "param.${i}.in" | tail -n5 | head -n1 | awk '{print $7,$8,$9,$10,$11,$12}' | sed 's|^|'${i}'\t|'
done
