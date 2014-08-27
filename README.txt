This program uses the InsightToolkit-4.6.0 with Module_SCIFIO=ON

To avoid SCIFIO throwhing exeption because of memory limits:
export JAVA_FLAGS=-Xmx3400m

Data to be analyzed for colocalization can be found here:
https://www.dropbox.com/sh/t2z0f4jttcnfh8s/AADq21HNR7EwH1JFxskhXXhea

Code for reading the files is taken from:
https://github.com/scifio/scifio-imageio/blob/master/test/itkSCIFIOImageIOTest.cxx

Use of itkSCIFIOImageIOTest with ome-tiff files:
./SCIFIOTestDriver itkSCIFIOImageIOTest bfconvert_output/output_series_0.ome.tiff test.ome.tiff -w -a -d 5

Known Bugs:
Mark Hiner (20.08.2014): "I'm assuming the # of planes to read isn't being updated for each series.. I'll look into it when I get a chance."
And I guess I have a memory leak (which might also be related to SCIFIO).
(https://github.com/scifio/scifio-imageio/issues)

Original .lif files have been converted with Bio-Formats 5.0.3 with:
./bfconvert /path/to/live-A.lif /path/to/live-A.ome.tiff

Bio-Formats 5.0.3 can be downloaded here:
http://downloads.openmicroscopy.org/bio-formats/5.0.3/


Some more code is taken from:
http://www.itk.org/Insight/Doxygen/html/IO_2ImageReadExtractWrite_8cxx-example.html
http://www.itk.org/Insight/Doxygen/html/Iterators_2ImageSliceIteratorWithIndex_8cxx-example.html
InsightToolkit/Examples/Statistics/ImageHistogram4.cxx











