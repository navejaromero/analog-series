#!/usr/bin/env bash

cat output/washed/* | tail -n +2 | cut -f2 > output/uwashed/tmp.txt
mkdir output/uwashed/TMP 
sort --temporary-directory=./output/uwashed/TMP --parallel=5 -u output/uwashed/tmp.txt > output/uwashed/uwashed.txt
rm output/uwashed/tmp.txt
rmdir output/uwashed/TMP
mkdir split
split -l 1000000 output/uwashed/uwashed.txt -d output/uwashed/split/uwashed_ --additional-suffix=".txt"
rm output/uwashed/uwashed.txt 
