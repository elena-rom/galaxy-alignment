# Galaxy Alignments

This is a series of programs used to determine galaxy clusters' principal axes from their moments of inertia, as well as comparing these to each cluster's brightest cluster galaxy's principal axis. 

It runs primarily on Python, using Astropy and Matplotlib. Data are used from the PanSTARRS-1 survey.

## Getting Started

Inside the working directory, create a directory titled "data" with all *.csv files containing galaxy coordinates for each cluster. Then, run script.sh to determine the moment of inertia for each file. The option -w will write all results to a file called result.csv, while -p will print all results to terminal.

Scrape.py can be used to download fits files centered at the BCG coordinate from the PanSTARR-1 survey. Several parameters must be set within this file.

Axes of BCGs should be manually recorded by running Source Extractor on each fits file. These should be recorded in a *.csv file, and then analysis.py can be run to compare the moment of inertia angles to the BCG angles.


## Prerequisites
* Python 3.x, Astropy, and Matplotlib for analysis

## Acknowledgments

Thanks to Michael West at Lowell Observatory and Roberto de Propris at the University of Turku for providing guidance and data for this project.
