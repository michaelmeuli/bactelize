#include <iostream>
#include <iomanip>
#include "itkImageToHistogramFilter.h"
#include "itkSCIFIOImageIO.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkRGBPixel.h"
#include "itkMetaDataObject.h"
#include "itkStreamingImageFilter.h"

int main( int argc, char * argv [] )
{

  if( argc < 2 )
    {
    std::cerr << "Missing command line arguments" << std::endl;
    std::cerr << "Usage :  Test input-ome-tiff-ImageFileName " << std::endl;
    return -1;
    }

  typedef unsigned short PixelComponentType;
  const unsigned int Dimension = 3;
  const unsigned int Channels = 3;
  typedef itk::Vector<PixelComponentType, Channels> PixelType;  
  typedef itk::Image<PixelType, Dimension> ImageType;

  typedef typename itk::ImageFileReader< ImageType > ReaderType;
  itk::SCIFIOImageIO::Pointer io = itk::SCIFIOImageIO::New();
  io->DebugOn();
  typename ReaderType::Pointer reader = ReaderType::New();
  std::cout << "reader->GetUseStreaming(): " << reader->GetUseStreaming() << std::endl;
  std::cout << "done checking streaming usage" << std::endl;
  reader->SetImageIO( io );
  const char * inputFileName  = argv[1];
  reader->SetFileName( inputFileName );

  typedef itk::StreamingImageFilter< ImageType, ImageType > StreamingFilter;
  typename StreamingFilter::Pointer streamer = StreamingFilter::New();
  streamer->SetInput( reader->GetOutput() );
  streamer->SetNumberOfStreamDivisions( 4 );


  try
    {
    streamer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Problem encoutered while reading image file : " << argv[1] << std::endl;
    std::cerr << excp << std::endl;
    return -1;
    }


  typedef itk::Statistics::ImageToHistogramFilter<ImageType>   HistogramFilterType;
  HistogramFilterType::Pointer histogramFilter = HistogramFilterType::New();
  typedef HistogramFilterType::HistogramSizeType   SizeType;
  SizeType size( 3 );

  size[0] =  40;        // number of bins for the Red   channel
  size[1] =   1;        // number of bins for the Green channel
  size[2] =   1;        // number of bins for the Blue  channel
  histogramFilter->SetHistogramSize( size );

  histogramFilter->SetMarginalScale( 10.0 ); 
  HistogramFilterType::HistogramMeasurementVectorType lowerBound( 3 );
  HistogramFilterType::HistogramMeasurementVectorType upperBound( 3 );
  lowerBound[0] = 0;
  lowerBound[1] = 0;
  lowerBound[2] = 0;
  upperBound[0] = 65536;
  upperBound[1] = 65536;
  upperBound[2] = 65536;
  histogramFilter->SetHistogramBinMinimum( lowerBound );
  histogramFilter->SetHistogramBinMaximum( upperBound ); 
  histogramFilter->SetInput(  streamer->GetOutput()  );

  histogramFilter->Update();
  
  typedef HistogramFilterType::HistogramType  HistogramType;
  const HistogramType * histogram = histogramFilter->GetOutput();
  const unsigned int histogramSize = histogram->Size();
  std::cout << std::endl << "Histogram size " << histogramSize << std::endl;
 
  unsigned int channel = 0;  // red channel
  std::cout << std::endl << "Histogram of the red component" << std::endl;
  for( unsigned int bin=0; bin < histogramSize; bin++ )
    {
    std::cout << "bin = " << std::setw(3) << bin << 
      "      measurement = " << std::setw(8) << std::setprecision(1) << std::setiosflags(std::ios::fixed) <<  histogram->GetMeasurement (bin, channel) <<
      "      frequency = " << std::setw(8) << histogram->GetFrequency( bin, channel ) << std::endl;	
    }

  size[0] =   1;  // number of bins for the Red   channel
  size[1] =  40;  // number of bins for the Green channel
  size[2] =   1;  // number of bins for the Blue  channel
  histogramFilter->SetHistogramSize( size );
  histogramFilter->Update();
  channel = 1;  // green channel
  std::cout << std::endl << "Histogram of the green component" << std::endl;
  for( unsigned int bin=0; bin < histogramSize; bin++ )
    {
    std::cout << "bin = " << std::setw(3) << bin << 
      "      measurement = " << std::setw(8) << std::setprecision(1) << std::setiosflags(std::ios::fixed) <<  histogram->GetMeasurement (bin, channel) <<
      "      frequency = " << std::setw(8) << histogram->GetFrequency( bin, channel ) << std::endl;	
    }

  size[0] =   1;  // number of bins for the Red   channel
  size[1] =   1;  // number of bins for the Green channel
  size[2] =  40;  // number of bins for the Blue  channel
  histogramFilter->SetHistogramSize( size );
  histogramFilter->Update();
  channel = 2;  // blue channel
  std::cout << std::endl << "Histogram of the blue component" << std::endl;
  for( unsigned int bin=0; bin < histogramSize; bin++ )
    {
    std::cout << "bin = " << std::setw(3) << bin << 
      "      measurement = " << std::setw(8) << std::setprecision(1) << std::setiosflags(std::ios::fixed) <<  histogram->GetMeasurement (bin, channel) <<
      "      frequency = " << std::setw(8) << histogram->GetFrequency( bin, channel ) << std::endl;	
    }
  std::cout << std::endl;

}
