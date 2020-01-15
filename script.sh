#!/bin/sh

# Instructions:
# If option -w is selected, will write results to file. If option -p is selected, will print results to terminal.
# Data csv files should be in a directory titled "data" in the working directory
# Will overwrite the result file each time it is run

while getopts ":wp" opt; do
  case $opt in
    w)
      echo "Writing results to file result.csv"
      > result.csv
      for filename in $PWD/data/*.csv; 
      do python3 momentofinertia.py ${filename} >> result.csv; 
      printf "\rProcessing $filename"
      done
      ;;
    p)
      for filename in $PWD/data/*.csv; do python3 momentofinertia.py ${filename}; done
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
  esac
done
