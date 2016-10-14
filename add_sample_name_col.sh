#!/bin/bash

shopt -s nullglob

for f in *.txt
do
  echo "Processing $f file..."
	sample=$(echo ${f} | cut -c1-5)
	awk -v sample=$sample '{print (sample, $0)}' < "$f" > "${f%.txt}.s.txt"
  echo "Finished $f file."
done
