#!/bin/sh

> result.csv
for filename in data/*.csv; do python3 momentofinertia2.py /home/elenarom/galaxy_alignment/"${filename}" >> result.csv; done