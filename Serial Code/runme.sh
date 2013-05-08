#!/bin/bash
FILES=../Images/*

for f in $FILES
do
	./surf $f ./Images/$f | tee -a data.csv
done
