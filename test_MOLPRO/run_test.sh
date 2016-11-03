#!/bin/bash

MOLP=molpros_2012.1
TMPDIR=$PWD

if [ ! -e mp2.out ]; then
 $MOLP -W$TMPDIR -d$TMPDIR -I$TMPDIR mp2
fi

### Extract equilibrium geometry in Newton-X format
awk 'BEGIN{take="none"}
     /Atomic Masses/{take="mass";n=NR;nat=0}
     /Atomic Coordinates/{take="geom"}
     /Frequencies dumped to record/{take="none"}
     {
      if (take =="mass" && NF==0) {take="none"}
      if (take =="mass" && NR>n) {for (i=3;i<=NF;i++){nat++;mass[nat]=$i}}
      if (take =="geom") {
       at[$1]=$2;noat[$1]=$3;x[$1]=$4;y[$1]=$5;z[$1]=$6
      }}
      END{
       for (i=1;i<=nat;i++){
        print at[i],noat[i],x[i],y[i],z[i],mass[i]
      }}' mp2.out  > geom

### Extract the Force Constants
awk '
     /Force Constants/{pr=1;n=NR}
     {if (NF==0) {pr=0}
      if (pr==1 && NR>n) {
      if ($2+0 != $2){for (i=1;i<=NF;i++){coord1[i]=$i;ncoord++;coord[ncoord]=$i}}
      else{for (i=2;i<=NF;i++){fc[$1,coord1[i-1]]=$i}}
     }}END{
     for (i=1;i<=ncoord;i++){
      for (j=1;j<=i;j++){
       printf " "fc[coord[i],coord[j]]
      } 
      printf "\n"
     }}' mp2.out > hessian.dat
       
../bin/sharc_init.exe < init.in > init.out
