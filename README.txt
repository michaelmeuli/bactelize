This program uses the InsightToolkit (compiled from most recent master on 2014/10/15) with Module_SCIFIO=ON

To avoid SCIFIO throwhing exeption because of memory limits:
export JAVA_FLAGS=-Xmx5400m

Data to be analyzed for colocalization can be found here:
https://www.dropbox.com/sh/t2z0f4jttcnfh8s/AADq21HNR7EwH1JFxskhXXhea

Code for reading the files is taken from:
https://github.com/scifio/scifio-imageio/blob/master/test/itkSCIFIOImageIOTest.cxx

Use of itkSCIFIOImageIOTest with ome-tiff files:
./SCIFIOTestDriver itkSCIFIOImageIOTest /path/to/dead-A.ome.tiff /path/to/dead-A.ome.tiff -w -a -d 5

Original .lif files have been converted with Bio-Formats 5.0.5 with:
./bfconvert /path/to/live-A.lif /path/to/live-A.ome.tiff

Bio-Formats 5.0.5 can be downloaded here:
http://downloads.openmicroscopy.org/bio-formats/5.0.5/











