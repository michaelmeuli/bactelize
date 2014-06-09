This program extensively uses the itk and vtk libraries: 
http://www.itk.org
http://www.vtk.org

How to build:
Build vtk firs and then itk with:
Module_ITKVtkGlue=ON
Module_SCIFIO=ON
(On a fresh Ubuntu 14.04 install you might want to do first:
sudo apt-get install subversion gcc g++ libX11-dev libXt-dev libgl1-mesa-dev libosmesa6-dev libglu1-mesa-dev git cmake cmake-curses-gui build-essential git-core)

To avoid SCIFIO throwhing exeption because of memory limits:
export JAVA_FLAGS=-Xmx3400m

Data to be analyzed for colocalization can be found here:
https://www.dropbox.com/sh/t2z0f4jttcnfh8s/AADq21HNR7EwH1JFxskhXXhea

Original .lif files have been converted with Bio-Formats 5.0.2 with:
./bfconvert /path/to/live-A.lif /path/to/live-A_%s.ome.tiff
http://downloads.openmicroscopy.org/bio-formats/5.0.2/

Lots of code is taken from examples found there:
InsightToolkit/Examples
http://www.itk.org/Wiki/ITK/Examples
https://github.com/scifio/scifio-imageio/tree/master/test
http://www.itk.org/Insight/Doxygen/html/IO_2ImageReadExtractWrite_8cxx-example.html
http://www.itk.org/Insight/Doxygen/html/Iterators_2ImageSliceIteratorWithIndex_8cxx-example.html
InsightToolkit/Examples/Statistics/ImageHistogram4.cxx



Use of itkSCIFIOImageIOTest with ome-tiff:
./SCIFIOTestDriver itkSCIFIOImageIOTest bfconvert_output/output_series_0.ome.tiff test.ome.tiff -w -d 5

SCIFIO was fixed Thu Apr 10 12:45:05 2014 -0500: 
commit 9f6f245a10b7837e6d7303b872ed55611fd99b5a
Author: Mark Hiner <hinerm@gmail.com>
Date:   Thu Apr 10 12:45:05 2014 -0500
    ENH: bump to latest scifio-imageio   
    Updated the SCIFIO-ImageIO to fix a pixel type detection error. It
    seemed that the returned pixel type was clashing with the
    UNKNOWNCOMPONENTTYPE constant, causing an unknown component error to be
    thrown erroneously.   
    Change-Id: I8f691f82da9288709e48626dd13d102d9df15bff







