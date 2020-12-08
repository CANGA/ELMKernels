#!/bin/bash 
awk '/mean/'  mpionly.results > trashfile
awk '{ print $3 $5 $7}' trashfile > trashfile1
awk '{ print $2}' resources.txt > trashfile2
paste -d "" trashfile2 trashfile1 > ../../../../results/utils/mpionly.dat

rm trashfile*
