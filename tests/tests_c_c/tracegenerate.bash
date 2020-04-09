#!/bin/bash 

awk -F 'mean:' '{ print $2 }' mpionly_smt1.log  >> dmpismt1.log
awk NF dmpismt1.log > dmpi_smt1.log

awk -F 'mean:' '{ print $2 }' mpionly_smt4.log  >> dmpismt4.log
awk NF dmpismt4.log > dmpi_smt4.log


rm dmpismt*.log

paste dmpi_smt1.log dmpi_smt4.log > data_mpionly.log
