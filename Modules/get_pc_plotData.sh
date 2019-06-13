#!/bin/bash

class=$1
dat=$2
odir=$3
sep=$4

string=$(sed -n -e "/^${dat} <- c/p" "${odir}${sep}${class}_model.r" | sed "s/^${dat} <- c(//g" | sed 's/)$//g')
vec=(${string//,/ })
printf '%s\n' "${vec[@]}" > "${odir}${sep}${class}_data_${dat}.txt"
