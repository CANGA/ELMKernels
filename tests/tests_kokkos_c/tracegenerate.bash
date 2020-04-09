#!/bin/bash 

awk -F 'mean:' '{ print $2 }' serialsmt1.log  >> dssmt1.log
awk NF dssmt1.log > dserial_smt1.log

awk -F 'mean:' '{ print $2 }' serialsmt4.log  >> dssmt4.log
awk NF dssmt4.log > dserial_smt4.log


rm dssmt*.log

paste dserial_smt1.log dserial_smt4.log > data_kokkos_serial.log
