#! /usr/bin/bash

MERGED_SRC_FILE="$1" 

OUTPUT_FILE="$2"

cdo -selvar,ivtnorm -expr,'ivtnorm=sqrt((1/9.81 * (husua^2)) + (1/9.81 * (husva^2)))' -vertsum -selvar,husua,husva -expr,'husua=hus*ua; husva=hus*va;' -ml2pl,100000,97500,95000,92500,90000,87500,85000,82500,80000,77500,75000,70000,65000,60000,55000,50000,45000,40000,35000,30000,25000,22500,20000,17500,15000,12500,10000,7000,5000,3000 $MERGED_SRC_FILE $OUTPUT_FILE

