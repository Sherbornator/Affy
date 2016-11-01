
#!/bin/bash

shopt -s nullglob

for f in *.vcf
do
	echo "Processing $f file..."
	python vcf_concatenate.py $f
	echo "Finished $f file."
done
