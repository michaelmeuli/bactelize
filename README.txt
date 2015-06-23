This program uses the "Insight Segmentation and Registration Toolkit (ITK) 4.7" with Module_SCIFIO=ON

To avoid SCIFIO throwhing exeption because of memory limits:
export JAVA_FLAGS=-Xmx5400m

Data to be analyzed for colocalization can be found here:
https://www.dropbox.com/sh/pe0aj43ts114b1p/AABhaM6KkUvczEkq1a73Op5Ia?dl=0

Known Bugs:
image2DReadMean is not always correct

Code for reading the files is taken from:
https://github.com/scifio/scifio-imageio/blob/master/test/itkSCIFIOImageIOTest.cxx

Original .lif files have been converted with Bio-Formats 5.0.5 with:
./bfconvert /path/to/live-A.lif /path/to/live-A.ome.tiff

Bio-Formats 5.0.5 can be downloaded here:
http://downloads.openmicroscopy.org/bio-formats/5.0.5/

Use of itkSCIFIOImageIOTest with ome-tiff files:
./SCIFIOTestDriver itkSCIFIOImageIOTest /path/to/dead-A.ome.tiff /path/to/dead-A.ome.tiff -w -a -d 5

On a fresh Ubuntu 14.04 install you might want to do first:
sudo apt-get install subversion gcc g++ libX11-dev libXt-dev libgl1-mesa-dev libosmesa6-dev libglu1-mesa-dev git cmake cmake-curses-gui build-essential git-core












