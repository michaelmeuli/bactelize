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
#include "itksys/SystemTools.hxx"
#include "QuickView.h"
#include "vnl/vnl_math.h"
#include "itkImageSliceConstIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"


int main( int argc, char * argv [] )
{

  if ( argc < 3 )
    {
    std::cerr << "Missing parameters. " << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
              << " inputImageFile outputImageFile"
              << std::endl;
    return -1;
    }




  typedef unsigned short              PixelType;
  typedef itk::Image< PixelType, 2 >  ImageType2D;
  typedef itk::Image< PixelType, 3 >  ImageType3D;

  typedef itk::ImageLinearIteratorWithIndex< ImageType2D > LinearIteratorType;
  typedef itk::ImageSliceConstIteratorWithIndex< ImageType3D > SliceIteratorType;

  typedef itk::ImageFileReader< ImageType3D > ReaderType;
  typedef itk::ImageFileWriter< ImageType2D > WriterType;


//  ReaderType::Pointer reader = ReaderType::New();
//  reader->SetFileName( argv[1] );


/*
  typedef unsigned short PixelComponentType;
  const unsigned int Dimension = 2;
  typedef unsigned int PixelType;
//  const unsigned int Channels = 3;
//  typedef itk::Vector<PixelComponentType, Channels> PixelType;  
  typedef itk::Image<PixelType, Dimension> ImageType;

  typedef typename itk::ImageFileReader< ImageType > ReaderType;
*/

  itk::SCIFIOImageIO::Pointer io = itk::SCIFIOImageIO::New();
  io->DebugOn();
  typename ReaderType::Pointer reader = ReaderType::New();
  std::cout << "reader->GetUseStreaming(): " << reader->GetUseStreaming() << std::endl;
  std::cout << "done checking streaming usage" << std::endl;
  reader->SetImageIO( io );
  const char * inputFileName  = argv[1];
  reader->SetFileName( inputFileName );

  typedef itk::StreamingImageFilter< ImageType3D, ImageType3D > StreamingFilter;
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


  ImageType3D::ConstPointer inputImage;
  inputImage = streamer->GetOutput();

//--
  unsigned int projectionDirection = 2;

  unsigned int i, j;
  unsigned int direction[2];
  for (i = 0, j = 0; i < 3; ++i )
    {
    if (i != projectionDirection)
      {
      direction[j] = i;
      j++;
      }
    }

  ImageType2D::RegionType region;
  ImageType2D::RegionType::SizeType size;
  ImageType2D::RegionType::IndexType index;

  ImageType3D::RegionType requestedRegion = inputImage -> GetRequestedRegion();

  index[ direction[0] ]    = requestedRegion.GetIndex()[ direction[0] ];
  index[ 1- direction[0] ] = requestedRegion.GetIndex()[ direction[1] ];
  size[ direction[0] ]     = requestedRegion.GetSize()[  direction[0] ];
  size[ 1- direction[0] ]  = requestedRegion.GetSize()[  direction[1] ];

  region.SetSize( size );
  region.SetIndex( index );

  ImageType2D::Pointer outputImage = ImageType2D::New();

  outputImage->SetRegions( region );
  outputImage->Allocate();
 
  SliceIteratorType  inputIt(  inputImage, inputImage->GetRequestedRegion() );
  LinearIteratorType outputIt( outputImage, outputImage->GetRequestedRegion() );

  inputIt.SetFirstDirection(  direction[1] );
  inputIt.SetSecondDirection( direction[0] );

  outputIt.SetDirection( 1 - direction[0] );

  outputIt.GoToBegin();
  while ( ! outputIt.IsAtEnd() )
    {
    while ( ! outputIt.IsAtEndOfLine() )
      {
      outputIt.Set( itk::NumericTraits<unsigned short>::NonpositiveMin() );
      ++outputIt;
      }
    outputIt.NextLine();
    }

  inputIt.GoToBegin();
  outputIt.GoToBegin();

  while( !inputIt.IsAtEnd() )
    {
    while ( !inputIt.IsAtEndOfSlice() )
      {
      while ( !inputIt.IsAtEndOfLine() )
        {
        outputIt.Set( vnl_math_max( outputIt.Get(), inputIt.Get() ));
        ++inputIt;
        ++outputIt;
        }
      outputIt.NextLine();
      inputIt.NextLine();

      }
    outputIt.GoToBegin();
    inputIt.NextSlice();
    }

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput(outputImage);
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
//--






  QuickView viewer;
  viewer.AddImage(outputImage.GetPointer(), true, itksys::SystemTools::GetFilenameName(argv[1]));  
  viewer.Visualize();



/* 

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
      "        measurement = " << std::setw(10) << std::setprecision(1) << std::setiosflags(std::ios::fixed) <<  histogram->GetMeasurement (bin, channel) <<
      "        frequency = " << std::setw(10) << histogram->GetFrequency( bin, channel ) << std::endl;	
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
      "        measurement = " << std::setw(10) << std::setprecision(1) << std::setiosflags(std::ios::fixed) <<  histogram->GetMeasurement (bin, channel) <<
      "        frequency = " << std::setw(10) << histogram->GetFrequency( bin, channel ) << std::endl;	
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
      "        measurement = " << std::setw(10) << std::setprecision(1) << std::setiosflags(std::ios::fixed) <<  histogram->GetMeasurement (bin, channel) <<
      "        frequency = " << std::setw(10) << histogram->GetFrequency( bin, channel ) << std::endl;	
    }
  std::cout << std::endl;

*/





}
