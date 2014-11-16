#! /bin/bash
  
set -e

for f in /home/michael/bioimage/batch/lif/*.lif; do
  echo "Processing " `basename $f` "..."
  /home/michael/bioimage/bftools-5.0.5/bfconvert "$f" /home/michael/bioimage/batch/ome/`basename "${f%.lif}.ome.tiff"`
done
