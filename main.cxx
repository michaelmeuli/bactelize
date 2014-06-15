/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/


#include "bactelize.h" 

int main( int argc, char * argv [] )
{
  if ( argc < 3 )
    {
    std::cerr << "Missing parameters. " << std::endl;
    std::cerr << "Usage: " << argv[0] << " series.ome.tiff-file" << " mip-projection-file-of-ch1" << std::endl;
    return -1;
    }

  int seriesnr = 0; 
  float xspacing = 0.33;
  float yspacing = 0.33;
  float zspacing = 1.2;
  float tspacing = 1.0;
  float cspacing = 1.0;
  int nucleichannel      = 0;
  int bacteriachannel    = 1;
  int lysotrackerchannel = 2;

  const char * inputFileName  = argv[1];
  SeriesReader seriesreader(inputFileName);
  std::cout << "Getting 5D Image of series number: " << seriesnr << std::endl;
  ImageType5D::Pointer image5D = seriesreader.get5DImage(seriesnr);
  seriesreader.dumpimageio();
  dumpmetadatadic(image5D);
  setspacing(image5D, xspacing, yspacing, zspacing, tspacing, cspacing);
  ImageType3D::ConstPointer inputImageMIP = extractchannel(image5D, bacteriachannel);


  // Histogram of second channell
  typedef itk::Statistics::ImageToHistogramFilter<ImageType3D>   HistogramFilterType;
  HistogramFilterType::Pointer histogramFilter = HistogramFilterType::New();
  typedef HistogramFilterType::HistogramSizeType   SizeType;
  SizeType size( 1 );
  size[0] =  40;        // number of bins for the green channel
  histogramFilter->SetHistogramSize( size );

  histogramFilter->SetMarginalScale( 10.0 ); 
  HistogramFilterType::HistogramMeasurementVectorType lowerBound( 3 );
  HistogramFilterType::HistogramMeasurementVectorType upperBound( 3 );
  lowerBound[0] = 0;
  upperBound[0] = 65536;
  histogramFilter->SetHistogramBinMinimum( lowerBound );
  histogramFilter->SetHistogramBinMaximum( upperBound ); 
  histogramFilter->SetInput(  inputImageMIP  );
  histogramFilter->Update();
  
  typedef HistogramFilterType::HistogramType  HistogramType;
  const HistogramType * histogram = histogramFilter->GetOutput();
  const unsigned int histogramSize = histogram->Size();
  std::cout << std::endl << "Histogram size " << histogramSize << std::endl;
 
  std::cout << std::endl << "Histogram of the green channell" << std::endl;
  for( unsigned int bin=0; bin < histogramSize; bin++ )
    {
    std::cout << "bin = " << std::setw(3) << bin << 
      "        measurement = " << std::setw(10) << std::setprecision(1) << std::setiosflags(std::ios::fixed) <<  histogram->GetMeasurement (bin, 0) <<
      "        frequency = " << std::setw(10) << histogram->GetFrequency( bin, 0 ) << std::endl;	
    }


  ImageType2D::Pointer outputImageMIP = maxintprojection(inputImageMIP);

  typedef itk::Image<unsigned char, 2>  ImageTypeWriter;
  typedef itk::ImageFileWriter< ImageTypeWriter > WriterType;
  typedef itk::RescaleIntensityImageFilter< ImageType2D, ImageTypeWriter >  RescaleFilterType;

  RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetInput( outputImageMIP ); 
  rescaleFilter->SetOutputMinimum( 0 );
  rescaleFilter->SetOutputMaximum( 255 );
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput(rescaleFilter->GetOutput());
  try
    {
    writer->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return -1;
    }




  QuickView viewer;
  viewer.AddImage(outputImageMIP.GetPointer(), true, itksys::SystemTools::GetFilenameName(argv[1]));  
  viewer.Visualize();

  return EXIT_SUCCESS;
}
