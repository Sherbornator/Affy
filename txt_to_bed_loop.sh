#!/bin/bash

shopt -s nullglob

for f in *.txt
do
  	echo "Processing $f file..."
        awk 'NR>1{printf("%s\t%s\t%s\t%s\n","chr"$2,$3,$3,$1":"$4":"$5)}' < "$f" > "${f%.txt}.bed"
        echo "Finished $f file."
done
