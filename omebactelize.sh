#! /bin/bash
  
set -e

for f in /home/michael/bioimage/batch/ome/*.ome.tiff; do
  echo "Processing " `basename $f` "..."
  /home/michael/bioimage/bactelize/build/bactelize /home/michael/bioimage/batch/ome/ /home/michael/bioimage/batch/out/ -f "$f" >> /home/michael/bioimage/batch/out/AA_term-out.txt 2>&1
done

exit 0




