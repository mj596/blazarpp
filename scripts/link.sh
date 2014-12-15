#! /bin/bash

path=`pwd`
for i in out/*; do
    echo ${path}/${i}
    rm -rf ${path}/${i}/LvPointAvg_QAccd 
    ln -s ${path}/${i}/LvTemplate_QAccd ${path}/${i}/LvPointAvg_QAccd 
done
