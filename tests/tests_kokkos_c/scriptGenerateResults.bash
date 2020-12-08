#!/bin/bash 

awk '/mean/' serial.results > trashfile
awk '{ print $3 $5 $7}' trashfile > trashfile1
awk '{ print $2}' resources.txt > trashfile2
paste -d "" trashfile2 trashfile1 > ../../../../results/utils/serial.dat

awk '/mean/' openmpsmt1.results > trashfile
awk '{ print $3 $5 $7}' trashfile > trashfile1
awk '{ print $2}' resources.txt > trashfile2
paste -d "" trashfile2 trashfile1 > ../../../../results/utils/openmpsmt1.dat

awk '/mean/' openmpsmt4.results > trashfile
awk '{ print $3 $5 $7}' trashfile > trashfile1
awk '{ print $2}' resources.txt > trashfile2
paste -d "" trashfile2 trashfile1 > ../../../../results/utils/openmpsmt4.dat

awk '/mean/' cuda.results > trashfile
awk '{ print $3 $5 $7}' trashfile > trashfile1
awk '{ print $1}' resources.txt | head -n 4 > trashfile2
paste -d "" trashfile2 trashfile1 > ../../../../results/utils/cuda.dat

rm trashfile*
