#! /bin/bash
 
set -e

/home/michael/bioimage/bactelize/build/bactelize /home/michael/bioimage/batch/ome/ /home/michael/bioimage/batch/out/ >> /home/michael/bioimage/batch/out/AA_logfile.txt 2>&1
