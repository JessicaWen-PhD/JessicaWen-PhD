#!/bin/bash
OUTPUT_NAME=Pt-Zundel-ontop
for i in {-10..10}; do
    cd ${OUTPUT_NAME}${i}
    sbatch run_scf.sh
    cd ../
done
